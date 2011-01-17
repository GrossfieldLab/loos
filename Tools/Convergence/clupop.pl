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
my $std_scale = 0;
my $zero_value = -20;


my $hdr = $0 . ' ' . join(' ', @ARGV);

my $ok = GetOptions(
		    'average!' => \$subtract_average,
		    'log!' => \$log_scale,
		    'std!' => \$std_scale,
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
  my $ravg = &columnAverage(\@M);
    
  for (my $j=0; $j<=$#M; ++$j) {
    for (my $i=0; $i<=$#bins; ++$i) {
      $M[$j]->[$i] -= $ravg->[$i];
    }
  }

  if ($std_scale) {
    my $rstd = &columnStd(\@M, $ravg);
    for (my $j=0; $j<=$#M; ++$j) {
      for (my $i=0; $i<=$#bins; ++$i) {
	$M[$j]->[$i] /= $rstd->[$i];
      }
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
Usage- clupop.pl [--average [--std]] [--log [--zero=value]] assignments >matrix.asc
EOF
  exit 0;
}


sub columnAverage {
  my $rM = shift;
  
  my @avgs = @{$rM->[0]};

  for (my $j=1; $j<=$#$rM; ++$j) {
    for (my $i=0; $i<=$#avgs; ++$i) {
      $avgs[$i] += $rM->[$j]->[$i];
    }
  }

  for (my $i=0; $i<=$#avgs; ++$i) {
    $avgs[$i] /= ($#$rM+1);
  }

  return(\@avgs);
}



sub columnStd {
  my $rM = shift;
  my $ravg = shift;

  my @stds = (0.0) x ($#$ravg+1);
  for (my $j=0; $j<=$#$rM; ++$j) {
    for (my $i=0; $i<=$#$ravg; ++$i) {
      my $d = $rM->[$j]->[$i] - $ravg->[$i];
      $stds[$i] += $d*$d;
    }
  }
  for (my $i=0; $i<=$#stds; ++$i) {
    $stds[$i] = sqrt($stds[$i]/$#$rM);
  }

  return(\@stds);
}
