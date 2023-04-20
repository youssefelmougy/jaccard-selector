#!/bin/bash

ROOT=$PWD/Selector
cd Selector

# Load modules
module load cray-openshmemx
module load cray-pmi
module unload darshan
module load perftools-base perftools

# Export environment variables
export PLATFORM=ex
export CC=cc
export CXX=CC
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/cray/pe/perftools/23.02.0/lib64

# Setup Bale
git clone https://github.com/jdevinney/bale.git bale
cd bale
patch -p1 < $ROOT/../scripts/perlmutter.patch
cd ..
cd bale/src/bale_classic
./bootstrap.sh
python3 ./make_bale -s
cd ../../../

# Setup HCLib
cd hclib
./install.sh --enable-production
source hclib-install/bin/hclib_setup_env.sh
cd modules/bale_actor && make
cd jaccard-selector
unzip ../inc/boost.zip -d ../inc/
export BALE_INSTALL=$ROOT/bale/src/bale_classic/build_${PLATFORM}
make
cd ../../../../
export BALE_INSTALL=$ROOT/bale/src/bale_classic/build_${PLATFORM}
export HCLIB_ROOT=$ROOT/hclib/hclib-install
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BALE_INSTALL/lib:$HCLIB_ROOT/lib:$HCLIB_ROOT/../modules/bale_actor/lib
export HCLIB_WORKERS=1


cd $ROOT/../

