#!/bin/bash

#SBATCH --job-name=compilazione
#SBATCH --partition=g100_usr_prod
#SBATCH --mem=3GB
#SBATCH --error  myJob.err
#SBATCH --output myJob.out
#SBATCH --account=pMI25_DCMIC
#SBATCH --time=00:10:00             # Tempo massimo stimato (HH:MM:SS)
#SBATCH --nodes=1                   # Numero di nodi
#SBATCH --ntasks-per-node=1         # Numero di task per nodo

# Carica il modulo di Singularity
module purge
module load profile/advanced
module load singularity

# compile.sh: Compiles the mf-adr project inside the Singularity container

CONTAINER="../hpc.sif"

# Prevent Singularity Stage 2 temp/cache errors by forcing it to scratch space
mkdir -p $CINECA_SCRATCH/singularity_tmp $CINECA_SCRATCH/singularity_cache
export SINGULARITY_TMPDIR=$CINECA_SCRATCH/singularity_tmp
export SINGULARITY_CACHEDIR=$CINECA_SCRATCH/singularity_cache

# Clean and prepare build directory
rm -rf build && mkdir build
cd build || exit 1

echo "Starting compilation inside $CONTAINER..."

# Run CMake and Make inside the container
singularity exec --cleanenv --bind $PWD/.. $CONTAINER bash -c "
    export PATH=/u/sw/toolchains/gcc-glibc/11.2.0/base/bin:$PATH
    export CC=mpicc
    export CXX=mpicxx
    cmake .. -DCMAKE_BUILD_TYPE=Release -DBoost_DIR="/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/trilinos/15.0.0/lib/external_packages/Boost" -Ddeal.II_DIR="/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/dealii/9.5.1/lib/cmake/deal.II" &&
    make -j \$(nproc)
"

echo "Compilation finished. Check ./build/single-node"
