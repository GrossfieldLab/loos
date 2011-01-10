#!/usr/bin/perl -w
#
# (c) 2011 Tod D. Romo, Grossfield Lab, URMC
#
# Cluster population over time (as a matrix) given a bin-assignment file...
#
#

use strict;

my @data;
my $maxbin = 0;
while (<>) {
  next if /^#/;
  chomp;
  push(@data, $_);
  if ($_ > $maxbin) {
    $maxbin = $_;
  }
}


++$maxbin;
my @bins = (0) x $maxbin;
my $n = 0;

for (my $i=0; $i<$#data; ++$i) {
  ++$bins[$data[$i]];
  ++$n;
  for (my $j=0; $j<=$#bins; ++$j) {
    print $bins[$j] / $n, "\t";
  }
  print "\n";
}
