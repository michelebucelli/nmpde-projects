#!/bin/bash

#SBATCH --job-name=compilazione_likwid          # Nome del job
#SBATCH --partition=g100_usr_prod
#SBATCH --mem=5GB
#SBATCH --error  myJob.err
#SBATCH --output myJob.out
#SBATCH --account=pMI25_DCMIC    # es: g100_all_serial, leonardo_usr_prod
#SBATCH --time=00:20:00             # Tempo massimo stimato (HH:MM:SS)
#SBATCH --nodes=1                   # Numero di nodi
#SBATCH --ntasks-per-node=1         # Numero di task per nodo

module purge
module load profile/advanced
module load singularity

# compile.sh: Compiles the mf-adr project inside the Singularity container

# Prevent Singularity Stage 2 temp/cache errors by forcing it to scratch space
export SINGULARITY_TMPDIR=$CINECA_SCRATCH/singularity_tmp
export SINGULARITY_CACHEDIR=$CINECA_SCRATCH/singularity_cache

mkdir -p $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR

# --- 1. DOWNLOAD AND BUILD LIKWID LOCALLY ---
export LIKWID_ROOT=$CINECA_SCRATCH/likwid_install
if [ ! -d "$LIKWID_ROOT" ]; then
    echo "perf_event backend..."
    
    cd $CINECA_SCRATCH/likwid_src
    
    # Patch config.mk for unprivileged HPC access
    sed -i "s|^PREFIX ?= .*|PREFIX = $LIKWID_ROOT|" config.mk
    sed -i "s|^ACCESSMODE = .*|ACCESSMODE = perf_event|" config.mk
    sed -i "s|^BUILDFREQ = .*|BUILDFREQ = false|" config.mk
    
    # Compile and install inside the container environment
    singularity exec --cleanenv -B $CINECA_SCRATCH /g100_scratch/userexternal/fghattas/hpc_proj/mf-adr/hpc.sif bash -c "
        source /u/sw/lmod/lmod/init/bash && 
	export MODULEPATH=/u/sw/toolchains/gcc-glibc/11.2.0/modules:/u/sw/modules/toolchains && 
	module load gcc-glibc dealii &&
        make && make install
    "
    cd -
fi
# ----------------------------------------------

rm -rf build && mkdir build
cd build || exit 1

# --- 2. COMPILE YOUR PROJECT WITH LIKWID ---
singularity exec --cleanenv -B $CINECA_SCRATCH /g100_scratch/userexternal/fghattas/hpc_proj/mf-adr/hpc.sif bash -c "
    export LD_LIBRARY_PATH=$LIKWID_ROOT/lib:\$LD_LIBRARY_PATH &&
    export PATH=/u/sw/toolchains/gcc-glibc/11.2.0/base/bin:$PATH
    export CC=mpicc
    export CXX=mpicxx
    cmake .. -DCMAKE_BUILD_TYPE=Release -DBoost_DIR="/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/trilinos/15.0.0/lib/external_packages/Boost" -Ddeal.II_DIR="/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/dealii/9.5.1/lib/cmake/deal.II" -DCMAKE_PREFIX_PATH=$LIKWID_ROOT -DCMAKE_CXX_FLAGS=\"-I$LIKWID_ROOT/include\" -DCMAKE_BUILD_RPATH=$LIKWID_ROOT/lib -DCMAKE_EXE_LINKER_FLAGS=\"-L$LIKWID_ROOT/lib -llikwid\" &&
    make -j \$(nproc)
"
