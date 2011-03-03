#!/usr/bin/perl -w
#
# Usage- bootstrap_overlap.pl model trajectory selection
#
# Notes:
# o If your BLAS/ATLAS are not in parallel, you can try running
#   concurrent BCOM/BOOTBCOM jobs with this script by using the
#   "--parallel" flag.
#
# Assumes the various LOOS tools are in your shell path
#
# Assumes you have gnuplot installed and is in your path
#
# The plotting is hard-coded.  If you want to change the plot, you
# will either need to change the code in the gnuplot function, or
# pattern your own gnuplot script after this.  Turning on debugging
# output (i.e. verbosity=3) will output the values used by the PERL
# program in gnuplot...
#
#
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



my $nreps = 10;        # Number of replicates for bootstrap
my $range = undef;     # Traj range to use (undef = auto)
my $stepsize = undef;  # Block size increment (undef = auto)
my $verbosity = 2;     # Extent of logging
my $prefix = undef;    # Intermediate file prefix (undef = auto)
my $scale = 1.0;       # Convert frame to time
my $units = 'ns';      # Units for plotting
my @seeds;             # Seeds for the exp fit in gnuplot
my $replot = 0;        # 1 = only generate plot
my $npts = 100;        # When range is auto, how many block sizes to use
my $parallel = 0;      # 1 = run boot/bcom concurrently

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
	       'parallel!' => \$parallel,
	       "help" => sub { &showHelp; }
	       
	      );

my $ok = GetOptions(%options);

&showHelp if (!$ok || $#ARGV != 2);


my $model = shift;
my $traj = shift;
my $sel = shift;


# Prefix is automatically generated from traj name, stripping off
# suffix
if (!defined($prefix)) {
  if ($traj =~ /\./) {
    $traj =~ /^(.+)\..+?$/;
    $prefix = $1;
  } else {
    $prefix = $traj;
  }
}



# Auto-generate the range if not specified by user.  We take the first
# half of the trajectory using the user-defined stepsize, or how many
# blocks were requested.  Stepsize has priority.

my $nframes = &getNumberOfFrames($model, $traj);
my $half = int($nframes/2);
if (!defined($range)) {
  my $delta = defined($stepsize) ? $stepsize : int($half / $npts);
  $range = "$delta:$delta:$half";
}


# If no seeds specified, use 1/2 and 1/20th of traj size
if ($#seeds < 0) {
  @seeds = ($half, $half/10.0);
}



if (!$replot) {


  if ($parallel) {
    # Run the BCOM job in the background if in parallel mode...
    my $child = &forkCommand("bcom --blocks $range $model $traj '$sel' >$prefix.bcom.asc");
    &runCommand("boot_bcom --blocks $range --replicates $nreps $model $traj '$sel' >$prefix.boot_bcom.asc");

    # Reap the forked proc...no zombies
    my $stat = waitpid $child, 0;

  } else {
    # Run jobs serially...
    &runCommand("bcom  --blocks $range $model $traj '$sel' >$prefix.bcom.asc");
    &runCommand("boot_bcom --blocks $range --replicates $nreps $model $traj '$sel' >$prefix.boot_bcom.asc");
  }
  
  # Combine the two separate output files from bcom and bootbcom into
  # one file for easy plotting...
  &mergeFiles("$prefix.conv.asc", "$prefix.bcom.asc", "$prefix.boot_bcom.asc");
}

# Call gnuplot to generate the plot
&gnuplot("$prefix.conv.ps", "$prefix.conv.asc", $prefix, $units, $scale, $nframes, \@seeds);



##############################################################################


# Forks, then runs command.
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
  my $fno = shift;    # output filename
  my $fni = shift;    # input filename
  my $ti = shift;     # base prefix name (used for labeling the plot)
  my $units = shift;  # Units string
  my $scale = shift;  # Scales frame into time
  my $num = shift;    # Size of traj
  my $rseeds = shift; # Seeds for exp fit

  # Get bounds for plot/fit
  my $half = int($num / 2);
  $num /= $scale;
  my($lower, $upper) = &findRange($fni);
  $lower /= $scale;
  $upper /= $scale;

  print "DEBUG> lower=$lower, upper=$upper\n" if ($verbosity > 2);
  print "DEBUG> half=$half, scale=$scale, units=$units\n" if ($verbosity > 2);
  print "DEBUG> seeds=(", join(',', @$rseeds), ")\n" if ($verbosity > 2);

  # Pipe output straight into gnuplot
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


# Run command providing logging (if necessary) and checks for errors
sub runCommand {
  my $cmd = shift;

  print STDERR "# $cmd\n" if ($verbosity > 1);
  my $failed = system($cmd);
  die "Error- '$cmd' failed with code ", $failed >> 8 if ($failed);
}


# Merges the BCOM & BOOTBCOM outputs into one file
# &mergeFiles(output, bcom, bootbcom)

sub mergeFiles {
  my $fno = shift;
  my $fn1 = shift;
  my $fn2 = shift;

  my $fho = new FileHandle "$fno", 'w'; defined($fho) || die "Error- cannot open $fno for output";
  my $fh1 = new FileHandle $fn1; defined($fh1) || die "Error- cannot open $fn1 for reading";
  my $fh2 = new FileHandle $fn2; defined($fh2) || die "Error- cannot open $fn2 for reading"; 

  my($f1, $f2);

  while (<$fh1>) {
    next if /^#/;     # Skip header metadata
    $f1 = $_;
    do {
      $f2 = <$fh2>;
    } while ($f2 =~ /^#/); # Also skip in the 2nd file
    chomp($f2);
    my @ary = split(/\s+/, $f2);
    shift(@ary);      # Shift off the frame no since it should be the
                      # same in both files
    chomp($f1);
    
    print $fho "$f1\t", join("\t", @ary), "\n";  # And merge
  }
}

# Finds the min and max time (frame no, in col-0) in the data file
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


# Calls trajinfo on a model/traj pair to get the # of frames in the
# trajectory
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


# Generate logging header
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
  my $seed_string = join(',', @seeds);
  my $parallel_string = $parallel ? 'yes' : 'no';
  my $range_string = defined($range) ? $range : 'auto';
  my $step_string = defined($stepsize) ? $stepsize : 'auto';
  my $prefix_string = defined($prefix) ? $prefix : 'auto';

  print <<EOF;
Usage- bootstrap_overlap.pl [options] model traj selection  
Options:
  --nreps=i ($nreps)\tNumber of replicates for bootstrap
  --range=s ($range_string)\tMatlab-style range of traj frames to use
  --stepsize=i ($step_string)\tBlock size increment
  --verbosity=i ($verbosity)\tAmount of logging
  --scale=f ($scale)\t\tScales frame index into time units
                    \t  e.g. units are in ns and each traj
                    \t  frame is in 0.1 ns, then scale=10
  --units=s ($units)\tString for time units in the plot
  --prefix=s ($prefix_string)\tPrefix for intermediate files
  --seed=f ($seed_string)\t\tSeeds for exponential fit in the plots
  --plotonly         \tOnly generate the plot (do not run BCOM/BOOTBCOM)
  --npts=i ($npts)   \tWhen auto-generating the traj range to use, use
                     \t  this number of points in the range (i.e. how
                     \t  many block-sizes to use
  --[no]parallel ($parallel_string)\tRun BCOM/BOOTBCOM concurrently
  --help
EOF

  exit(0);
}
