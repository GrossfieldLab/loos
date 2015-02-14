#!/bin/bash
#
# (c) 2015 Tod D. Romo, Grossfield Lab, URMC
#
# Build a tarball for MacOS suitable for distribution...
#
#

MAKE="make -j4"
SCONS="scons -j4"
CHECK=""

PREFIX=${1:-/Volumes/Doubler}
PACKAGES=${2:-/Volumes/Doubler/Packages}
CLEAN=${3:-}

TOPDIR=`pwd`

BOOSTDIR=`ls -lUd $PACKAGES/boost* | head -1 | awk '{print $9}'`
BOOST=`basename $BOOSTDIR`
echo "Using BOOST $BOOST"

NETCDFDIR=`ls -lUd $PACKAGES/netcdf* | head -1 | awk '{print $9}'`
NETCDF=`basename $NETCDFDIR`
echo "Using NetCDF $NETCDF"

HDF5DIR=`ls -lUd $PACKAGES/hdf5* | head -1 | awk '{print $9}'`
HDF5=`basename $HDF5DIR`
echo "Using HDF5 $HDF5"




# Determine LOOS version & install location
RELEASE_VERSION=`grep loos_version loos_build_config.py | cut -d\' -f2`
INSTALL_DIR="$PREFIX/loos-$RELEASE_VERSION"

if [ -d $INSTALL_DIR ] ; then
    echo "***WARNING***"
    echo "Install location already exists"
    echo "Moving it out of the way..."
    rm -f $INSTALL_DIR.bak
    mv $INSTALL_DIR $INSTALL_DIR.bak
fi
mkdir $INSTALL_DIR


# First, handle boost
echo '***MARKER: BOOST'
pushd $BOOSTDIR
# Really need a distclean here...
./bootstrap.sh --prefix=$INSTALL_DIR/$BOOST || exit -1
./b2 install
popd

# Now, deal with HDF5
echo '***MARKER: HDF5'
pushd $HDF5DIR
if [ -n "$CLEAN" ] ; then make distclean  ; fi
./configure --prefix=$INSTALL_DIR/$HDF5 || exit -1
$MAKE $CHECK || exit -1
$MAKE install || exit -1
popd

# Now, handle netcdf
echo '***MARKER: NETCDF'
pushd $NETCDFDIR
if [ -n "$CLEAN" ] ; then make distclean  ; fi
CPPFLAGS="-I$INSTALL_DIR/$HDF5/include" LDFLAGS="-L$INSTALL_DIR/$HDF5/lib" ./configure --prefix=$INSTALL_DIR/$NETCDF || exit -1
$MAKE $CHECK install || exit -1
popd

echo '***MARKER: LOOS'
# Now, build LOOS
if [ -n "$CLEAN" ] ; then scons -c ; scons -c config ; fi
$SCONS PREFIX=$INSTALL_DIR BOOST=$INSTALL_DIR/$BOOST NETCDF=$INSTALL_DIR/$NETCDF docs || exit -1
$SCONS PREFIX=$INSTALL_DIR BOOST=$INSTALL_DIR/$BOOST NETCDF=$INSTALL_DIR/$NETCDF install || exit -1

cd $PREFIX
tar cvf - loos-$RELEASE_VERSION | bzip2 -cv9 >$HOME/loos-$RELEASE_VERSION-macos.tar.bz2
