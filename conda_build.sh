#!/bin/bash
#  This file is part of LOOS.
#
#  LOOS (Lightweight Object-Oriented Structure library)
#  Copyright (c) 2013, Tod D. Romo
#  Department of Biochemistry and Biophysics
#  School of Medicine & Dentistry, University of Rochester
#
#  This package (LOOS) is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation under version 3 of the License.
#
#  This package is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#  Alan Grossfield, University of Rochester, 12/6/2019


# USAGE: ./conda_build.sh ENVNAME
#   will create conda environment ENVNAME if it doesn't exist, or install
#   packages into ENVNAME if it already exists, then run scons to build
#   LOOS.  By default, it runs scons with 4 compile processes, but if you
#   want to change that number (e.g. to 8 or 1), you can supply that as a
#   second argument, e.g.
#   ./conda_build.sh ENVNAME 8

numprocs=4
envname="loos"
while getopts "yhj:ie:" opt; do
    case ${opt} in
        e )
            envname=$OPTARG
            echo "Will use conda env $envname"
            ;;
        j )
            numprocs=$OPTARG
            echo "Number of processors: $num"
            ;;
        i )
            echo "Will install LOOS in the conda environment"
            do_install=1
            ;;
        y )
            echo "Will install non-interactively"
            non_interactive=1
            ;;
        h )
            echo "Usage:"
            echo "    -h         Display this message"
            echo "    -i         Install LOOS into the conda env"
            echo "    -e NAME    Use conda env NAME"
            echo "    -j N       Use N processors while compiling"
            echo "    -y         Install non-interactively"
            exit 0
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            exit 0
            ;;
    esac
done

platform=`uname`

echo "Setting channel priority to strict"
conda config --set channel_priority strict

packages="python=3 swig scons numpy scipy boost openblas libnetcdf lapack compilers eigen"

# does this env exist (there must be a smarter way to do this)
envs=$(conda env list | awk '{print $1}' )
for i in $envs; do
    if [ "$i" = "$envname" ]; then
        found=1
        break
    fi
done

# Build up the conda installation command line
if [ -z $found ]; then
    echo "Creating conda environment $envname"
    command="conda create "
    #conda create -n $envname -c conda-forge $packages
else
    echo "Installing into existing environment $envname"
    command="conda install "
    #conda install -n $envname -c conda-forge $packages
fi

if [ -z $non_interactive ]; then
    command+="--yes "
fi

command+="-n $envname -c conda-forge $packages "

# Run the conda command
echo $command
eval $command

echo "Activating environment $envname"
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate $envname

if [ "$platform" = "Linux" ]; then
    export CXX=`which g++`
elif [ "$platform" = 'Darwin' ]; then
    export CXX=`which clang++`
else
    echo "Unknown platform $platform, assuming g++"
    export CXX=`which g++`
fi
echo "CXX set to $CXX"

if [ ${do_install} ]; then
    scons -j$numprocs PREFIX=$CONDA_PREFIX
    scons PREFIX=$CONDA_PREFIX -j1 install
else
    scons -j$numprocs
fi
