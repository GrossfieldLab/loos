#!/usr/bin/perl -w
#
# (c) 2011 Tod D. Romo, Grossfield Lab, URMC
#
# Cluster population over time (as a matrix) given a bin-assignment file...
#
#

use strict;
use Getopt::Long;


my $subtract_average = 0;
my $log_scale = 0;
my $zero_value = -20;


my $hdr = $0 . ' ' . join(' ', @ARGV);

my $ok = GetOptions(
		    'average!' => \$subtract_average,
		    'log!' => \$log_scale,
		    'zero=f' => \$zero_value,
		    'help' => sub { &help; }
		   );

$ok || die "Unknown option";

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
my @M;

for (my $i=0; $i<$#data; ++$i) {
  ++$bins[$data[$i]];
  ++$n;
  my @row;
  for (my $j=0; $j<=$#bins; ++$j) {
    push(@row, $bins[$j] / $n);
  }
  push(@M, \@row);
}

if ($subtract_average) {
  for (my $i=0; $i<=$#bins; ++$i) {
    my $avg = 0.0;
    for (my $j=0; $j<=$#M; ++$j) {
      $avg += $M[$j]->[$i];
    }
    $avg /= ($#M + 1);
    
    for (my $j=0; $j<=$#M; ++$j) {
      $M[$j]->[$i] -= $avg;
    }
  }
}

if ($log_scale) { 
  for (my $j=0; $j<=$#M; ++$j) {
    for (my $i=0; $i<=$#bins; ++$i) {
      if ($M[$j]->[$i] == 0) {
	$M[$j]->[$i] = $zero_value;
      } else {
	$M[$j]->[$i] = log(abs($M[$j]->[$i]));
      }
    }
  }
}


print "# $hdr\n";
for (my $j=0; $j<=$#M; ++$j) {
  for (my $i=0; $i<=$#bins; ++$i) {
    print $M[$j]->[$i], "\t";
  }
  print "\n";
}



sub help {
  print <<EOF;
Usage- clupop.pl [--average] [--log [--zero=value]] assignments >matrix.asc
EOF
  exit 0;
}
