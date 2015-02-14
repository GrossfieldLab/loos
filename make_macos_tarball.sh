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
PACKAGES="/Volumes/Doubler/Packages"

PREFIX=${1:-/Volumes/Doubler}
CLEAN="$2"



TOPDIR=`pwd`

# Getting the basename is primarily in case we want to install into package-specific dirs
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
    rm -rf $INSTALL_DIR.bak
    mv $INSTALL_DIR $INSTALL_DIR.bak
fi
mkdir $INSTALL_DIR


# First, handle boost
echo '***MARKER: BOOST'
pushd $BOOSTDIR
# Really need a distclean here...
BOOSTOPTS=""
if [ -n "$CLEAN" ] ; then BOOSTOPTS="-a" ; fi
./bootstrap.sh --prefix=$INSTALL_DIR || exit -1
./b2 $BOOSTOPTS install
popd

# Now, deal with HDF5
echo '***MARKER: HDF5'
pushd $HDF5DIR
if [ -n "$CLEAN" ] ; then make distclean  ; fi
./configure --prefix=$INSTALL_DIR || exit -1
$MAKE $CHECK || exit -1
$MAKE install || exit -1
popd

# Now, handle netcdf
echo '***MARKER: NETCDF'
pushd $NETCDFDIR
if [ -n "$CLEAN" ] ; then make distclean  ; fi
CPPFLAGS="-I$INSTALL_DIR/include" LDFLAGS="-L$INSTALL_DIR/lib" ./configure --prefix=$INSTALL_DIR || exit -1
$MAKE $CHECK install || exit -1
popd

echo '***MARKER: LOOS'
# Now, build LOOS
if [ -n "$CLEAN" ] ; then scons -c ; scons -c config ; fi
$SCONS PREFIX=$INSTALL_DIR BOOST=$INSTALL_DIR NETCDF=$INSTALL_DIR docs || exit -1
$SCONS PREFIX=$INSTALL_DIR BOOST=$INSTALL_DIR NETCDF=$INSTALL_DIR install || exit -1

cd $PREFIX
tar cvf - loos-$RELEASE_VERSION | bzip2 -cv9 >$HOME/loos-$RELEASE_VERSION-macos.tar.bz2
