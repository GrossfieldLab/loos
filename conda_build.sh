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


# USAGE:   ./conda_build.sh -e loos -j 8 -i
#          This will create a conda env called loos (or install in it if it
#          exists). It will run scons with 8 jobs, and will do an install into
#          the conda env when it's done.
#          ./conda_build -h
#          will show all options.

numprocs=4
envname="loos"
conda_install_command="conda"

while getopts "yhj:ie:c:" opt; do
    case ${opt} in
        e )
            envname=$OPTARG
            echo "Will use conda env: $envname"
            ;;
        j )
            numprocs=$OPTARG
            echo "Number of processors: $num"
            ;;
        i )
            echo "Will install LOOS in the current conda environment"
            do_install=1
            ;;
        y )
            echo "Will install non-interactively"
            non_interactive=1
            ;;
        c)
            conda_install_command=$OPTARG
            echo "Conda install/create command: ${conda_install_command}"
            ;;
        h )
            echo "Usage:"
            echo "    -h         Display this message"
            echo "    -i         Install LOOS into the conda env"
            echo "    -e NAME    Use conda env NAME"
            echo "    -j N       Use N processors while compiling"
            echo "    -y         Install non-interactively"
            echo "    -c NAME   Use name (e.g. mamba) for creating/installing"
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

packages="python=3 swig=4 cmake numpy scipy scikit-learn boost openblas libnetcdf lapack compilers eigen "

env_found=$(conda env list | egrep -v '^#' | egrep "^${envname}[ ]" )
# Build up the conda installation command line
if [[ ${env_found} ]] ; then
    echo "Installing into existing environment $envname"
    command="${conda_install_command} install "
else
    echo "Creating conda environment $envname"
    command="${conda_install_command} create "
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

if [[ -d build ]] ; then
    echo "Warning- using existing build directory"
else
    mkdir build
fi

cd build
echo "*** Configuring LOOS ***"
cmake -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} ..

echo "*** Building LOOS ***"
cmake --build . -j${numprocs}

if [[ ${do_install} ]]; then
    echo  "*** Installing LOOS ***"
    cmake --install .
fi
