#ifndef PRECONDITIONERS_HPP
#define PRECONDITIONERS_HPP

#include "./IncludeFiles.hpp"

using namespace dealii;

// ==================================================================
// Preconditioners for the unsteady Navier-Stokes saddle point system
//
// Block System Structure:
//     [  F   B^T ] [ U ]   [ G ]
//     [ -B    0  ] [ P ] = [ 0 ]
//
//  - F   : Momentum block (Velocity stiffness/mass/convection)
//  - B^T : Discrete gradient operator
//  - -B  : Discrete divergence operator
//  - U   : Velocity vector
//  - P   : Pressure vector
// ==================================================================

enum class Preconditioner
{
    IDENTITY,
    BLOCK_IDENTITY,
    SIMPLE,
    APPROX_SIMPLE,
    YOSIDA,
    APPROX_YOSIDA,
    PCD,
    APPROX_PCD
};

// ---------------------------------------------------------------
// Class: PreconditionIdentity
// Identity preconditioner for scalar Trilinos vectors: P^-1 v = v.
// ---------------------------------------------------------------
class PreconditionIdentity
{
public:
    void vmult(TrilinosWrappers::MPI::Vector &dst,
               const TrilinosWrappers::MPI::Vector &src) const
    {
        dst = src;
    }
};

// ---------------------------------------------------------------
// Class: PreconditionBlockIdentity
// Identity preconditioner for block Trilinos vectors: P^-1 v = v.
// ---------------------------------------------------------------
class PreconditionBlockIdentity
{
public:
    void vmult(TrilinosWrappers::MPI::BlockVector &dst,
               const TrilinosWrappers::MPI::BlockVector &src) const
    {
        dst = src;
    }
};

// ---------------------------------------------------------------
// Abstract Base Class: BlockPrecondition
// Base interface for block preconditioners with helper routines
// for inner solver initialization (ILU or AMG).
// ---------------------------------------------------------------
class BlockPrecondition
{
public:
    virtual ~BlockPrecondition() = default;

    /**
     * @brief Apply the block preconditioner: dst = P^{-1} * src.
     */
    virtual void vmult(TrilinosWrappers::MPI::BlockVector &dst,
                       const TrilinosWrappers::MPI::BlockVector &src) const = 0;

protected:
    /**
     * @brief Helper method to set up inner preconditioners (ILU or AMG)
     *        for individual sub-blocks.
     * @param preconditioner Shared pointer target for the initialized preconditioner.
     * @param matrix System matrix block to precondition.
     * @param ilu True to use Trilinos ILU, false for Trilinos AMG.
     */
    void initialize_inner_preconditioner(
        std::shared_ptr<TrilinosWrappers::PreconditionBase> &preconditioner,
        const TrilinosWrappers::SparseMatrix &matrix,
        bool ilu)
    {
        if (ilu)
        {
            auto actual_preconditioner =
                std::make_shared<TrilinosWrappers::PreconditionILU>();
            actual_preconditioner->initialize(matrix);
            preconditioner = actual_preconditioner;
        }
        else
        {
            auto actual_preconditioner =
                std::make_shared<TrilinosWrappers::PreconditionAMG>();
            actual_preconditioner->initialize(matrix);
            preconditioner = actual_preconditioner;
        }
    }
};

// ---------------------------------------------------------------
// Class: PreconditionSIMPLE
//
// Semi-Implicit Method for Pressure-Linked Equations (SIMPLE).
// Approximates the Schur complement S using the diagonal of F:
//     S_SIMPLE = -B * diag(F)^-1 * B^T
//
// Sign conventions expected:
//     negB_ : -B  (Divergence block scaled by -1)
//     B_t_  :  B^T (Gradient block)
// ---------------------------------------------------------------
class PreconditionSIMPLE : public BlockPrecondition
{
public:
    /**
     * @brief Initializes matrix pointers, computes negDinv = -diag(F)^{-1},
     *        assembles S = negB * diag(F)^{-1} * B_t, and builds inner preconditioners.
     */
    void initialize(const TrilinosWrappers::SparseMatrix &F_,
                    const TrilinosWrappers::SparseMatrix &negB_,
                    const TrilinosWrappers::SparseMatrix &B_t_,
                    const TrilinosWrappers::MPI::BlockVector &vec,
                    const double &alpha_ = 0.5,
                    const unsigned int &maxit_ = 10000,
                    const double &tol_ = 1e-2,
                    const bool &ilu = true)
    {
        F = &F_;
        negB = &negB_;
        B_t = &B_t_;

        alpha = alpha_; // Relaxation parameter in (0, 1]
        maxit = maxit_;
        tol = tol_;

        // Compute negDinv = -1 / diag(F)
        negDinv.reinit(vec.block(0));
        for (unsigned int index : negDinv.locally_owned_elements())
        {
            negDinv[index] = -1.0 / F->diag_element(index);
        }

        // Assemble approximate Schur complement matrix:
        // S = negB * (-negDinv) * B_t = -B * diag(F)^-1 * B^T
        negB->mmult(S, *B_t, negDinv);

        // Build inner preconditioners for F and S
        this->initialize_inner_preconditioner(preconditioner_F, *F, ilu);
        this->initialize_inner_preconditioner(preconditioner_S, S, ilu);
    }

    /**
     * @brief Exact SIMPLE application using inner GMRES solves for F and S.
     */
    void vmult(TrilinosWrappers::MPI::BlockVector &dst,
               const TrilinosWrappers::MPI::BlockVector &src) const override
    {
        tmp.reinit(src);

        // Step 1: Solve velocity block F * sol1_u = src_u
        SolverControl solver_control_F(maxit, tol * src.block(0).l2_norm());
        SolverGMRES<TrilinosWrappers::MPI::Vector> solver_F(solver_control_F);
        solver_F.solve(*F, tmp.block(0), src.block(0), *preconditioner_F);

        // Step 2: Solve pressure block S * sol1_p = B * sol1_u - src_p
        B_t->Tvmult(tmp.block(1), tmp.block(0));
        tmp.block(1) -= src.block(1);

        SolverControl solver_control_S(maxit, tol * tmp.block(1).l2_norm());
        SolverGMRES<TrilinosWrappers::MPI::Vector> solver_S(solver_control_S);
        solver_S.solve(S, dst.block(1), tmp.block(1), *preconditioner_S);

        // Step 3: Apply SIMPLE pressure relaxation
        dst.block(1) /= alpha;

        // Step 4: Correct velocity: dst_u = sol1_u - D^-1 * B^T * dst_p
        dst.block(0) = tmp.block(0);
        B_t->vmult(tmp.block(0), dst.block(1));
        tmp.block(0).scale(negDinv);
        dst.block(0) += tmp.block(0);
    }

protected:
    double alpha;
    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *negB;
    const TrilinosWrappers::SparseMatrix *B_t;
    TrilinosWrappers::MPI::Vector negDinv;
    TrilinosWrappers::SparseMatrix S;
    std::shared_ptr<TrilinosWrappers::PreconditionBase> preconditioner_F;
    std::shared_ptr<TrilinosWrappers::PreconditionBase> preconditioner_S;
    mutable TrilinosWrappers::MPI::BlockVector tmp;
    unsigned int maxit;
    double tol;
};

// ---------------------------------------------------------------
// Class: PreconditionApproxSIMPLE
// Inexact SIMPLE variant using single preconditioner applications
// (ILU/AMG) in place of inner GMRES iterations for speed.
// ---------------------------------------------------------------
class PreconditionApproxSIMPLE : public PreconditionSIMPLE
{
public:
    void vmult(TrilinosWrappers::MPI::BlockVector &dst,
               const TrilinosWrappers::MPI::BlockVector &src) const override
    {
        tmp.reinit(src);

        // Inexact solve: single application of preconditioner_F
        preconditioner_F->vmult(tmp.block(0), src.block(0));

        // Inexact solve: single application of preconditioner_S
        B_t->Tvmult(tmp.block(1), tmp.block(0));
        tmp.block(1) -= src.block(1);
        preconditioner_S->vmult(dst.block(1), tmp.block(1));

        // Pressure relaxation
        dst.block(1) /= alpha;

        // Velocity correction step
        dst.block(0) = tmp.block(0);
        B_t->vmult(tmp.block(0), dst.block(1));
        tmp.block(0).scale(negDinv);
        dst.block(0) += tmp.block(0);
    }
};

using PreconditionaSIMPLE = PreconditionApproxSIMPLE;

// ---------------------------------------------------------------
// Class: PreconditionYosida
//
// Yosida Algebraic Splitting Preconditioner.
// Approximates the Schur complement using the diagonal velocity mass matrix M_u:
//     S_Yosida = -dt * B * diag(M_u)^-1 * B^T
//
// Sign conventions expected:
//     negB_ : -B  (Divergence block scaled by -1)
//     B_t_  :  B^T (Gradient block)
// ---------------------------------------------------------------
class PreconditionYosida : public BlockPrecondition
{
public:
    /**
     * @brief Initializes Yosida preconditioner matrix blocks and computes
     *        Dinv = dt / diag(M_u).
     * @param dt_ Physical time step size (delta t).
     */
    void initialize(const TrilinosWrappers::SparseMatrix &F_,
                    const TrilinosWrappers::SparseMatrix &negB_,
                    const TrilinosWrappers::SparseMatrix &B_t_,
                    const TrilinosWrappers::SparseMatrix &M_u_,
                    const TrilinosWrappers::MPI::BlockVector &vec,
                    const double &dt_ = 1.0,
                    const unsigned int &maxit_ = 10000,
                    const double &tol_ = 1e-2,
                    const bool &ilu = true)
    {
        F = &F_;
        negB = &negB_;
        B_t = &B_t_;

        maxit = maxit_;
        tol = tol_;
        dt = dt_;

        // Compute Dinv = dt / diag(M_u)
        Dinv.reinit(vec.block(0));
        for (unsigned int index : Dinv.locally_owned_elements())
        {
            Dinv[index] = dt / M_u_.diag_element(index);
        }

        // Assemble Schur complement matrix: negS = -B * (dt * M_u^-1) * B^T
        negB->mmult(negS, *B_t, Dinv);

        // Build inner preconditioners for F and negS
        this->initialize_inner_preconditioner(preconditioner_F, *F, ilu);
        this->initialize_inner_preconditioner(preconditioner_S, negS, ilu);
    }

    /**
     * @brief Exact Yosida application with inner GMRES iterations.
     */
    void vmult(TrilinosWrappers::MPI::BlockVector &dst,
               const TrilinosWrappers::MPI::BlockVector &src) const override
    {
        tmp1.reinit(src);

        // Step 1: Solve F * sol1_u = src_u
        SolverControl solver_control_F(maxit, tol * src.block(0).l2_norm());
        SolverGMRES<TrilinosWrappers::MPI::Vector> solver_F(solver_control_F);
        solver_F.solve(*F, tmp1.block(0), src.block(0), *preconditioner_F);

        // Step 2: Solve -S * sol1_p = src_p - B * sol1_u
        tmp1.block(1) = src.block(1);
        negB->vmult_add(tmp1.block(1), tmp1.block(0));

        SolverControl solver_control_S(maxit, tol * tmp1.block(1).l2_norm());
        SolverGMRES<TrilinosWrappers::MPI::Vector> solver_S(solver_control_S);
        solver_S.solve(negS, dst.block(1), tmp1.block(1), *preconditioner_S);

        // Step 3: Second velocity solve for Yosida correction
        dst.block(0) = tmp1.block(0);
        B_t->vmult(tmp1.block(0), dst.block(1));

        tmp2.reinit(src.block(0));
        SolverControl solver_control_F2(maxit, tol * tmp1.block(0).l2_norm());
        SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_F2(solver_control_F2);
        solver_gmres_F2.solve(*F, tmp2, tmp1.block(0), *preconditioner_F);

        dst.block(0) -= tmp2;
    }

protected:
    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *negB;
    const TrilinosWrappers::SparseMatrix *B_t;
    TrilinosWrappers::MPI::Vector Dinv;
    TrilinosWrappers::SparseMatrix negS;
    std::shared_ptr<TrilinosWrappers::PreconditionBase> preconditioner_F;
    std::shared_ptr<TrilinosWrappers::PreconditionBase> preconditioner_S;
    mutable TrilinosWrappers::MPI::BlockVector tmp1;
    mutable TrilinosWrappers::MPI::Vector tmp2;
    unsigned int maxit;
    double dt;
    double tol;
};

// ---------------------------------------------------------------
// Class: PreconditionApproxYosida
// Inexact Yosida variant using single preconditioner applications
// (ILU/AMG) in place of inner GMRES solves.
// ---------------------------------------------------------------
class PreconditionApproxYosida : public PreconditionYosida
{
public:
    void vmult(TrilinosWrappers::MPI::BlockVector &dst,
               const TrilinosWrappers::MPI::BlockVector &src) const override
    {
        tmp1.reinit(src);

        // Inexact step 1: velocity block
        preconditioner_F->vmult(tmp1.block(0), src.block(0));

        // Inexact step 2: pressure block
        tmp1.block(1) = src.block(1);
        negB->vmult_add(tmp1.block(1), tmp1.block(0));
        preconditioner_S->vmult(dst.block(1), tmp1.block(1));

        // Inexact step 3: velocity correction
        dst.block(0) = tmp1.block(0);
        B_t->vmult(tmp1.block(0), dst.block(1));

        tmp2.reinit(src.block(0));
        preconditioner_F->vmult(tmp2, tmp1.block(0));

        dst.block(0) -= tmp2;
    }
};

using PreconditionaYosida = PreconditionApproxYosida;

// ---------------------------------------------------------------
// Class: PreconditionPCD
//
// Pressure Convection-Diffusion (PCD) Preconditioner.
// Approximates the inverse Schur complement as:
//     S_PCD^-1 = - A_p^-1 * F_p * M_p^-1
//
// Blocks required:
//  - Ap : Pressure Laplacian operator
//  - Mp : Pressure Mass matrix
//  - Fp : Pressure Convection-Diffusion operator
// ---------------------------------------------------------------
class PreconditionPCD : public BlockPrecondition
{
public:
    void initialize(const TrilinosWrappers::SparseMatrix &F_,
                    const TrilinosWrappers::SparseMatrix &B_t_,
                    const TrilinosWrappers::SparseMatrix &Ap_,
                    const TrilinosWrappers::SparseMatrix &Mp_,
                    const TrilinosWrappers::SparseMatrix &Fp_,
                    const unsigned int &maxit_ = 10000,
                    const double &tol_ = 1e-2,
                    const bool &ilu = true)
    {
        F = &F_;
        B_t = &B_t_;
        Ap = &Ap_;
        Mp = &Mp_;
        Fp = &Fp_;

        maxit = maxit_;
        tol = tol_;

        // Build inner preconditioners for F, Ap, and Mp
        this->initialize_inner_preconditioner(preconditioner_F, *F, ilu);
        this->initialize_inner_preconditioner(preconditioner_Ap, *Ap, ilu);
        this->initialize_inner_preconditioner(preconditioner_Mp, *Mp, ilu);
    }

    /**
     * @brief Exact PCD application using inner GMRES solves for Mp, Ap, and F.
     */
    void vmult(TrilinosWrappers::MPI::BlockVector &dst,
               const TrilinosWrappers::MPI::BlockVector &src) const override
    {
        tmp.reinit(src);

        // Step 1: q1 = Mp^-1 * src_p
        TrilinosWrappers::MPI::Vector q1(src.block(1));
        SolverControl solver_control_Mp(maxit, tol * src.block(1).l2_norm());
        SolverGMRES<TrilinosWrappers::MPI::Vector> solver_Mp(solver_control_Mp);
        solver_Mp.solve(*Mp, q1, src.block(1), *preconditioner_Mp);

        // Step 2: q2 = Fp * q1
        TrilinosWrappers::MPI::Vector q2(src.block(1));
        Fp->vmult(q2, q1);

        // Step 3: q3 = Ap^-1 * q2
        TrilinosWrappers::MPI::Vector q3(src.block(1));
        SolverControl solver_control_Ap(maxit, tol * q2.l2_norm());
        SolverGMRES<TrilinosWrappers::MPI::Vector> solver_Ap(solver_control_Ap);
        solver_Ap.solve(*Ap, q3, q2, *preconditioner_Ap);

        dst.block(1) = q3;
        dst.block(1) *= -1.0;

        // Step 4: Solve for velocity block F * dst_u = src_u - B^T * dst_p
        TrilinosWrappers::MPI::Vector Btp(src.block(0));
        B_t->vmult(Btp, dst.block(1));
        tmp.block(0) = src.block(0);
        tmp.block(0) -= Btp;

        SolverControl solver_control_F(maxit, tol * tmp.block(0).l2_norm());
        SolverGMRES<TrilinosWrappers::MPI::Vector> solver_F(solver_control_F);
        solver_F.solve(*F, dst.block(0), tmp.block(0), *preconditioner_F);
    }

protected:
    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *B_t;
    const TrilinosWrappers::SparseMatrix *Ap;
    const TrilinosWrappers::SparseMatrix *Mp;
    const TrilinosWrappers::SparseMatrix *Fp;
    std::shared_ptr<TrilinosWrappers::PreconditionBase> preconditioner_F;
    std::shared_ptr<TrilinosWrappers::PreconditionBase> preconditioner_Ap;
    std::shared_ptr<TrilinosWrappers::PreconditionBase> preconditioner_Mp;
    mutable TrilinosWrappers::MPI::BlockVector tmp;
    unsigned int maxit;
    double tol;
};

// ---------------------------------------------------------------
// Class: PreconditionApproxPCD
// Inexact PCD variant using single preconditioner applications
// (ILU/AMG) for Mp, Ap, and F.
// ---------------------------------------------------------------
class PreconditionApproxPCD : public PreconditionPCD
{
public:
    void vmult(TrilinosWrappers::MPI::BlockVector &dst,
               const TrilinosWrappers::MPI::BlockVector &src) const override
    {
        tmp.reinit(src);

        // Inexact Mp solve
        TrilinosWrappers::MPI::Vector q1(src.block(1));
        preconditioner_Mp->vmult(q1, src.block(1));

        // Pressure convection application
        TrilinosWrappers::MPI::Vector q2(src.block(1));
        Fp->vmult(q2, q1);

        // Inexact Ap solve
        TrilinosWrappers::MPI::Vector q3(src.block(1));
        preconditioner_Ap->vmult(q3, q2);

        dst.block(1) = q3;
        dst.block(1) *= -1.0;

        // Inexact velocity solve
        TrilinosWrappers::MPI::Vector Btp(src.block(0));
        B_t->vmult(Btp, dst.block(1));
        tmp.block(0) = src.block(0);
        tmp.block(0) -= Btp;

        preconditioner_F->vmult(dst.block(0), tmp.block(0));
    }
};

using PreconditionaPCD = PreconditionApproxPCD;

// ---------------------------------------------------------------
// Class: PreconditionBlockTriangular
//
// Generic Upper or Lower Block Triangular Preconditioner given a
// custom user-provided Schur complement approximation matrix.
// ---------------------------------------------------------------
/*
class PreconditionBlockTriangular : public BlockPrecondition
{
public:
    void initialize(const TrilinosWrappers::SparseMatrix &F_,
                    const TrilinosWrappers::SparseMatrix &negB_,
                    const TrilinosWrappers::SparseMatrix &B_t_,
                    const TrilinosWrappers::SparseMatrix &Schur_approx_,
                    const unsigned int &maxit_ = 10000,
                    const double &tol_ = 1e-2,
                    const bool &ilu = true,
                    const bool &is_upper_ = true)
    {
        F = &F_;
        neg_B = &negB_;
        B_t = &B_t_;
        Schur_approx = &Schur_approx_;

        maxit = maxit_;
        tol = tol_;
        is_upper = is_upper_;

        this->initialize_inner_preconditioner(preconditioner_F, *F, ilu);
        this->initialize_inner_preconditioner(preconditioner_S, *Schur_approx, ilu);
    }

    void vmult(TrilinosWrappers::MPI::BlockVector &dst,
               const TrilinosWrappers::MPI::BlockVector &src) const override
    {
        tmp.reinit(src);

        if (is_upper)
        {
            // Solve Schur block first, then velocity block
            SolverControl solver_control_S(maxit, tol * src.block(1).l2_norm());
            SolverGMRES<TrilinosWrappers::MPI::Vector> solver_S(solver_control_S);
            solver_S.solve(*Schur_approx, dst.block(1), src.block(1), *preconditioner_S);

            TrilinosWrappers::MPI::Vector Btp(src.block(0));
            B_t->vmult(Btp, dst.block(1));
            tmp.block(0) = src.block(0);
            tmp.block(0) -= Btp;

            SolverControl solver_control_F(maxit, tol * tmp.block(0).l2_norm());
            SolverGMRES<TrilinosWrappers::MPI::Vector> solver_F(solver_control_F);
            solver_F.solve(*F, dst.block(0), tmp.block(0), *preconditioner_F);
        }
        else
        {
            // Solve velocity block first, then Schur block
            SolverControl solver_control_F(maxit, tol * src.block(0).l2_norm());
            SolverGMRES<TrilinosWrappers::MPI::Vector> solver_F(solver_control_F);
            solver_F.solve(*F, dst.block(0), src.block(0), *preconditioner_F);

            tmp.block(1) = src.block(1);
            neg_B->vmult_add(tmp.block(1), dst.block(0));

            SolverControl solver_control_S(maxit, tol * tmp.block(1).l2_norm());
            SolverGMRES<TrilinosWrappers::MPI::Vector> solver_S(solver_control_S);
            solver_S.solve(*Schur_approx, dst.block(1), tmp.block(1), *preconditioner_S);
        }
    }

private:
    bool is_upper;
    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *neg_B;
    const TrilinosWrappers::SparseMatrix *B_t;
    const TrilinosWrappers::SparseMatrix *Schur_approx;
    std::shared_ptr<TrilinosWrappers::PreconditionBase> preconditioner_F;
    std::shared_ptr<TrilinosWrappers::PreconditionBase> preconditioner_S;
    mutable TrilinosWrappers::MPI::BlockVector tmp;
    unsigned int maxit;
    double tol;
};
*/

#endif // PRECONDITIONERS_HPP