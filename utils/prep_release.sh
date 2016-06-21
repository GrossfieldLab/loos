#!/bin/bash
#
# (c) 2008 Tod D. Romo
#     Grossfield Lab
#     URMC
#
#  NOTE:  This script is used internally to create a release from the
#         Git repository



export PREF=/var/tmp/loos-inst   # Where to put test install of LOOS
export RELDIR=/var/tmp           # Where to export LOOS release for tarballing
export PROCS=4                   # Number of procs to use in scons build

export GIT=git                   # Change this to echo for debugging (disables checkins)

######################################

function AbortOrContinue {
    local OK

    echo $1
    read -p 'Do you want to proceed?  [y/N] ' OK
    if [ "$OK" = "N" -o "$OK" = "n" -o \( -z "$OK" \) ] ; then
	echo "Release aborted!"
	exit -1
    fi
}


######################################

GITDIR=`pwd`

echo "*** Checking GIT status"
if ! $GIT diff-index --quiet HEAD -- ; then
    echo "ERROR - you appear to have unchecked modifications to this release."
    echo "        Please fix this and run this script again..."
    echo
    $GIT diff --name-status
    exit -1
fi




echo "Checking version ID..."
DOXVERS=`grep '^PROJECT_NUMBER' Doxyfile | gawk '{print $3}' | sed 's/^v//'`
SCONSVERS=`egrep '^loos_version' loos_build_config.py | gawk '{print $3}' | sed "s/'//g"`
if [ "$SCONSVERS" != "$DOXVERS" ] ; then
    echo "Error- scons version is $SCONSVERS, but Doxyfile version is $DOXVERS"
    exit -1
fi

AbortOrContinue "You are building release for version $SCONSVERS ?"



VERS=$SCONSVERS

pushd $RELDIR
$GIT clone $GITDIR loos-$VERS
cd loos-$VERS
rm -rf .git

echo "*** Building documentation"
echo "+ Cleaning..."
rm -rf Doc                # Manually remove since scons wont
scons -cs caboodle        # Clean everything (to be safe)
echo "+ Building..."
scons -sj$PROCS docs      # Rebuild docs explicitly

echo "*** Checking install target"
scons -sj$PROCS install PREFIX=$PREF
CWD=`pwd`
if ./utils/check_loos_install.pl --nofull --exclude '.git' --exclude 'utils' --exclude 'membrane_map.hpp' --prep `pwd` $PREF ; then
    echo "Install appears ok.  Generating full report..."
    ( echo "***INSTALL APPEARS OK***" ;\
      echo "Please check the diff list below for any errors..." ;\
      echo ;\
      ./utils/check_loos_install.pl --full `pwd` $PREF ) | less

    AbortOrContinue

else
    echo "***SUSPECTED ERRORS IN INSTALLATION***"
    echo "Aborting"
    exit -1
fi

echo "*** Removing test install"
rm -r $PREF




cd ..
tar cvf - loos-$VERS | gzip -cv9 >~/loos-$VERS.tar.gz
rm -r loos-$VERS
