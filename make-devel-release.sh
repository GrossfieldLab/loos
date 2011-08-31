#!/bin/bash -x
#
# (c) 2011 Tod D. Romo, Grossfield Lab, URMC
#
#

# Divine release number
REVISION=`svn info | perl -ne 'chomp;if(/^URL/){@a=split(m@/@,$_);$b=$a[$#a];}if(/^Revision: (\d+)/){$r=$1;}END{print"loos-$b-$r\n";}'`

# Sanitize the changelog
cp ChangeLog ChangeLog-
sed -r 's/@[^>]+//' <ChangeLog- >ChangeLog
rm ChangeLog-


tar cvf - --transform="s/^\./$REVISION/" --exclude=.svn --exclude='*~' --exclude='.scon*' . | bzip2 -cv9 >$HOME/$REVISION.tar.bz2

svn revert ChangeLog

