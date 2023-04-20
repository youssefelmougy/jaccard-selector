#!/bin/bash

# Load modules
module load cray-openshmemx
module load cray-pmi
module unload darshan
module load perftools-base perftools

export PLATFORM=ex
export BALE_INSTALL=$PWD/Selector/bale/src/bale_classic/build_${PLATFORM}
export HCLIB_ROOT=$PWD/Selector/hclib/hclib-install
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BALE_INSTALL/lib:$HCLIB_ROOT/lib:$HCLIB_ROOT/../modules/bale_actor/lib
export HCLIB_WORKERS=1
export CC=cc
export CXX=CC
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/cray/pe/perftools/23.02.0/lib64

module load cmake/3.22.0
module load gcc/10.3.0

