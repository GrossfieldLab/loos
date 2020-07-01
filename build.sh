#!/bin/sh
# This install script is intended for conda-forge, and assumes the conda env
# is already set up. If you need to set up a conda env, I suggest running
# conda_build.sh instead.


echo "CXX set to $CXX"

scons PREFIX=$CONDA_PREFIX
scons PREFIX=$CONDA_PREFIX install
