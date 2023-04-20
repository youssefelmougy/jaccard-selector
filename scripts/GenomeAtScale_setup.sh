#!/bin/bash

ROOT=$PWD/GenomeAtScale
cd GenomeAtScale

# Load modules
module load cmake/3.22.0
module load gcc/10.3.0

# Setup BLAS library (v3.10.0)
if [ ! -d "BLAS-3.10.0" ];
then
    wget http://www.netlib.org/blas/blas-3.10.0.tgz
    tar xvfz blas-3.10.0.tgz
    cd BLAS-3.10.0
    cp make.inc make2.inc
    sed 's/-frecursive/-frecursive -fPIC/g' make2.inc >! make.inc
    make -j
    ln -s blas_LINUX.a libblas.a
    cd ..
fi

# Setup lapack library (v3.10.1)
if [ ! -d "lapack-3.10.1" ];
then
    wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.1.tar.gz
    tar xvfz v3.10.1.tar.gz
    cd lapack-3.10.1
    sed 's/-O\([0-9]\)\+/-O\1 -fPIC/g' make.inc.example | sed 's/LDFLAGS =/LDFLAGS = -fPIC/g' > make.inc
    make -j
    cd ..
fi

# Setup CTF library
cd ctf
FFLAGS="-fallow-argument-mismatch" ./configure LIBS='-lblas -llapack -lgfortran' LIB_PATH="-L$ROOT/BLAS-3.10.0 -L$ROOT/lapack-3.10.1" --no-dynamic CXX="CC -DNULL=0" --with-lapack --with-scalapack
make -j
cd scalapack
sed 's/mpif90/ftn -fallow-argument-mismatch/g' SLmake.inc.example | sed 's/mpicc/cc/g' > SLmake.inc
make BLASLIB="-L$ROOT/BLAS-3.10.0 -lblas" LAPACKLIB="-L$ROOT/lapack-3.10.1 -llapack"

# Setup jaccard-ctf code (GenomeAtScale)
cd ../../jaccard-ctf
make clean
make CXX=CC INCLUDE_PATH="-I$ROOT/ctf/include" LIB_PATH="-L$ROOT/ctf/lib -L$ROOT/BLAS-3.10.0 -L$ROOT/lapack-3.10.1 -L$ROOT/ctf/scalapack" LIB_FILES="-lctf -lscalapack -llapack -lblas"
unzip kmer_matrix.zip


cd $ROOT/../

