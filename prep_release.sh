#!/bin/bash
#
# (c) Tod D. Romo
#     Grossfield Lab
#     URMC
#
#  NOTE:  This file is NOT part of any release!
#
#

### Small sanity check...
OK=`pwd | egrep 'tags/release-[0-9.]+$'`
if [ ! $OK ] ; then
    echo "ERROR - you don't appear to be prepping a release in a release directory!"
    echo "        Make sure you are in a properly tagged release from the SVN"
    echo "        prior to running this script."
    exit -1
fi


# Build the documentation...

scons -cs caboodle
rm -rf Docs

# One more sanity check...
OK=`svn -q status`
if [ -n "$OK" ] ; then
    echo "ERROR - you appear to have unchecked modifications to this release."
    echo "        Please fix this and run this script again..."
    echo
    echo "> svn -q status"
    echo "$OK"
    exit -1
fi

scons -s docs

# Massage the ChangeLog

cp ChangeLog ChangeLog-
sed -r 's/@[^>]+//' <ChangeLog- >ChangeLog
rm ChangeLog-


# Update SVN

RELEASE=`pwd | sed -r 's@^/.+/@@'`

read -p 'Shall I automatically update the SVN for you? (this will also remove this script) [Y/n] ' OK

if [ "$OK" = "Y" -o "$OK" = "y" -o \( -z "$OK" \) ] ; then
    echo "Updating SVN for $RELEASE..."
    
    set -o xtrace
    svn add Docs
    svn ci -m 'Added Documentation to $RELEASE' Docs
    svn ci -m 'Sanitized ChangeLog prior to $RELEASE' ChangeLog
    svn rm prep_release.sh
    svn ci -m 'Removed prep_release script from $RELEASE' prep_release.sh

else
    echo "***WARNING***WARNING***WARNING***"
    echo
    echo "DON'T FORGET YOU MUST UPDATE THE SVN PRIOR TO RELEASING $RELEASE"
fi


