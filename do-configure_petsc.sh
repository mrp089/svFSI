#!/bin/bash

build_dir=build
rm -rf $build_dir
mkdir $build_dir
cd $build_dir

# built from source
export PETSC_DIR=/Users/pfaller/work/repos/petsc/arch-darwin-c-opt

export LD_LIBRARY_PATH=$PETSC_DIR/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$PETSC_DIR/lib:$LIBRARY_PATH

cmake \
-DSV_USE_PETSC=ON \
.. 

make -j12 CC=gcc-13 FC=gfortran-13 CXX=g++-13
