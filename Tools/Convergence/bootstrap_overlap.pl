#!/usr/bin/perl -w
#
# Usage- bootstrap_overlap.pl model trajectory selection
#


#  This file is part of LOOS.
#
#  LOOS (Lightweight Object-Oriented Structure library)
#  Copyright (c) 2010, Tod D. Romo
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



use FileHandle;
use Getopt::Long;



my $nreps = 10;
my $range = undef;
my $stepsize = undef;
my $verbosity = 2;
my $prefix = undef;
my $scale = 1.0;
my $units = 'ns';
my @seeds;
my $replot = 0;
my $npts = 100;
my $parallel = 1;

my $hdr = &header($0, \@ARGV);


my %options = (
	       'nreps=i' => \$nreps,
	       'range=s' => \$range,
	       'stepsize=i' => \$stepsize,
	       'verbosity=i' => \$verbosity,
	       'scale=f' => \$scale,
	       'units=s' => \$units,
	       'prefix=s' => \$prefix,
	       'seed=f' => \@seeds,
	       'plotonly' => \$replot,
	       'npts=i' => \$npts,
	       'parallel!' => \$parallel
	      );

my $ok = GetOptions(%options);

&showHelp(\%options) if (!$ok || $#ARGV != 2);


my $model = shift;
my $traj = shift;
my $sel = shift;


if (!defined($prefix)) {
  if ($traj =~ /\./) {
    $traj =~ /^(.+)\..+?$/;
    $prefix = $1;
  } else {
    $prefix = $traj;
  }
}


my $nframes = &getNumberOfFrames($model, $traj);
my $half = int($nframes/2);
if (!defined($range)) {
  my $delta = defined($stepsize) ? $stepsize : int($half / $npts);
  $range = "$delta:$delta:$half";
}

if ($#seeds < 0) {
  @seeds = ($half, $half/10.0);
}



if (!$replot) {
  if ($parallel) {
    my $child = &forkCommand("bcom $model $traj '$sel' 1 $range >$prefix.bcom.asc");
    &runCommand("boot_bcom $model $traj '$sel' $nreps 0 1 $range >$prefix.boot_bcom.asc");
    my $stat = waitpid $child, 0;
  } else {
    &runCommand("bcom $model $traj '$sel' 1 $range >$prefix.bcom.asc");
    &runCommand("boot_bcom $model $traj '$sel' $nreps 0 1 $range >$prefix.boot_bcom.asc");
  }
  
  &mergeFiles("$prefix.conv.asc", "$prefix.bcom.asc", "$prefix.boot_bcom.asc");
}
&gnuplot("$prefix.conv.ps", "$prefix.conv.asc", $prefix, $units, $scale, $nframes, \@seeds);



##############################################################################


sub forkCommand {
  my $cmd = shift;

  print "# $cmd\n" if ($verbosity > 1);
  my $pid = fork;
  
  if (!defined($pid)) {
    die "Error- cannot fork off '$cmd'";
  } elsif ($pid == 0) {
    my $ok = system("$cmd");
    exit($ok);
  } else {
    return($pid);
  }

}



sub gnuplot {
  my $fno = shift;
  my $fni = shift;
  my $ti = shift;
  my $units = shift;
  my $scale = shift;
  my $num = shift;
  my $rseeds = shift;

  my $half = int($num / 2);
  $num /= $scale;
  my($lower, $upper) = &findRange($fni);
  $lower /= $scale;
  $upper /= $scale;

  print STDERR "DEBUG> lower=$lower, upper=$upper\n";

  my $fh = new FileHandle "|gnuplot";
  defined($fh) || die 'Error- cannot open pipe to gnuplot';
  print $fh <<EOF;
set out "$fno"
set term post portrait enhanced color

set title "$ti - $num $units"
set xlabel "Block size ($units)"

set style line 1 lc 1 lw 1 lt 1
set style line 2 lc 1 lw 2 lt 1

set style line 3 lc 2 lw 1 lt 1
set style line 4 lc 2 lw 2 lt 1

set style line 5 lc 3 lw 2 lt 1
set style line 6 lc 4 lw 4 lt 2

a=0.5
b=$$rseeds[0]
c=0.5
d=$$rseeds[1]


f(x) = a*exp(-x/b)+c*exp(-x/d)+e*exp(-x/f)+1
a=1;b=10;c=1;d=100;e=1;f=1000
fit f(x) "$fni" u (\$1/$scale):(\$5/\$2) via a,b,c,d,e,f

set multiplot layout 2,1
set ylabel "Covariance Overlap"
set yrange [0.3:1]
set xrange [$lower:$upper]

set key bottom right

plot "$fni" every 2::1 u (\$1/$scale):2:(sqrt(\$3)) w error ls 1 notitle,\\
"" u (\$1/$scale):2:(sqrt(\$3)) w l ls 2 ti 'BCOM',\\
\\
"" every 2::1 u (\$1/$scale):5:(sqrt(\$6)) w error ls 3 notitle,\\
"" u (\$1/$scale):5:(sqrt(\$6)) w l ls 4 title '(R=$nreps) BOOT BCOM'


unset title
set ylabel "Bootstrap / Blocked"
set yrange [*:*]
set key top right

plot "$fni" u (\$1/$scale):(\$5/\$2) w l ls 5 ti 'Ratio',\\
"" u (\$1/$scale):(f(\$1/$scale)) w l ls 6 ti sprintf('y=%.3g exp(-t/%.3g) + %.3g exp(-t/%.3g) + %.3g exp(-t/%.3g) + 1', a, b, c, d, e, f)

unset multiplot

set out "error_$fno"
set term post enhanced solid color landscape
set xrange [*:*]
set yrange [*:*]
g(x) = 0
plot "$fni" u (\$1/$scale):(\$5/\$2 - f(\$1/$scale)) w l ti 'Residual', g(x) w l not
quit

EOF
}


sub runCommand {
  my $cmd = shift;

  print STDERR "# $cmd\n" if ($verbosity > 1);
  my $failed = system($cmd);
  die "Error- '$cmd' failed with code ", $failed >> 8 if ($failed);
}


sub mergeFiles {
  my $fno = shift;
  my $fn1 = shift;
  my $fn2 = shift;

  my $fho = new FileHandle "$fno", 'w'; defined($fho) || die "Error- cannot open $fno for output";
  my $fh1 = new FileHandle $fn1; defined($fh1) || die "Error- cannot open $fn1 for reading";
  my $fh2 = new FileHandle $fn2; defined($fh2) || die "Error- cannot open $fn2 for reading"; 

  my($f1, $f2);

  while (<$fh1>) {
    next if /^#/;
    $f1 = $_;
    do {
      $f2 = <$fh2>;
    } while ($f2 =~ /^#/);
    chomp($f2);
    my @ary = split(/\s+/, $f2);
    shift(@ary);
    chomp($f1);
    
    print $fho "$f1\t", join("\t", @ary), "\n";
  }
}


sub findRange {
  my $fn = shift;
  my $fh = new FileHandle $fn; defined($fh) || die "Error- cannot open $fn";
  
  my $mint = 1e100;
  my $maxt = 0.0;

  while (<$fh>) {
    next if /^#/;
    chomp;
    my @ary = split;
    if ($ary[0] < $mint) {
      $mint = $ary[0];
    }
    if ($ary[0] > $maxt) {
      $maxt = $ary[0];
    }
  }

  return($mint, $maxt);
}


sub getNumberOfFrames {
  my $model = shift;
  my $traj = shift;

  my $fh = new FileHandle "trajinfo -b $model $traj 2>&1|";
  defined($fh) || die "Error- cannot open pipe from trajinfo command";

  my $dummy = <$fh>;
  my @ary = split(/\s+/, $dummy);
  $#ary == 3 || die "Error- trajinfo command gave bad results";
  return($ary[1]);
}


sub header {
  my $prog = shift;
  my $ra = shift;

  my $hdr = $prog . ' ';
  foreach (@$ra) {
    $hdr .= '\'' . $_ . '\' ';
  }
  chop($hdr);
  return($hdr);
}




sub showHelp {
  my $rh = shift;

  print "Usage- bootstrap_overlap.pl [options] model traj selection\n";
  print "Options:\n";
  foreach (sort keys %$rh) {
    next if ($_ eq 'help');

    my $rval = $rh->{$_};
    my $val;

    # Is this option bound to a variable?
    if (ref($rval) ne '') {
      # Special handling of arrays
      if (ref($rval) eq 'ARRAY') {
	$val = '(' .  join(',', @$rval) . ')';
      } else {
	if (!defined($$rval)) {
	  $val = '<not set>';
	} else {
	  $val = "'$$rval'";
	}
      }
    } else {
      die "Should not be here!";
    }

    printf "\t%20s = %s\n", $_, $val;
  }


  exit(0);
}
