# Setup environment for bourne-shell-like shells...
#
#  This file is part of LOOS.
#
#  LOOS (Lightweight Object-Oriented Structure library)
#  Copyright (c) 2012 Tod D. Romo, Grossfield Lab
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

###
### Set this to where LOOS is installed, or where the uninstalled distribution is...
###
DEFAULT_PATH=""
###
###
###

CWD=`pwd`
LOOSPATH="${DEFAULT_PATH:-$CWD}"


# Find library
if [ -d $LOOSPATH/lib ] ; then
    LIBPATH="$LOOSPATH/lib"
else
    LIBPATH="$LOOSPATH"
fi

# Find tools
if [ -d $LOOSPATH/bin ] ; then
    TOOLPATH="$LOOSPATH/bin"
else
    TOOLPATH="$LOOSPATH/Tools:$LOOSPATH/Packages/Convergence:$LOOSPATH/Packages/DensityTools:$LOOSPATH/Packages/ElasticNetworks:$LOOSPATH/Packages/HydrogenBonds:$LOOSPATH/Packages/User"
fi

export PATH="$PATH:$TOOLPATH"

# Use DYLD_LIBRARY_PATH or LD_LIBRARY_PATH depending on OS
SYSTEM=`uname`
if [ "$SYSTEM" = "Darwin" ] ; then
    export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$LIBPATH"
else
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$LIBPATH"
fi

