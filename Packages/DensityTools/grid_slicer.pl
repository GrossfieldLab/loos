#!/usr/bin/perl -w
#
# (c) 2009 Tod D. Romo
#          Grossfield Lab
#          University of Rochester Medical and Dental School
#
# Extracts slices from a grid and makes figures of them...
#
#


use FileHandle;
use Getopt::Long;

$stride = 0;
$indices = '';
$help_flag = 0;
$gridbase = '';

$result = GetOptions("stride=i" => \$stride,
		     "range=s" => $indices,
		     "output=s" => \$gridbase,
		     "help" => \$help_flag
		    );

$result || die "$0: Error parsing command line options...";
!($stride && ($indices ne '')) || die "$0: Error - you cannot specify both a stride and a range";
if (!$stride && $indices eq '') {
  $stride = 10;
}
&show_help if ($help_flag || $#ARGV < 0);

$gridname = shift;
if ($gridbase eq '') {
  $gridbase = $gridname; $gridbase =~ s/\.[^.]*$//;
}

@dims = &getGridSize($gridname);
print "Grid size is (", join(',', @dims), ")\n";

if ($stride) {
  $indices = "0:$stride:$dims[2]";
}
@indices = &parseRanges($indices);
print "Selecting ", $#indices+1, " slices.\n";

print STDERR 'Processing- ';
$tilecount = 0;
@gnuscript = ( 'set term post enhanced color solid' );
foreach $index (@indices) {
  my $basename = sprintf "%s_%03d", $gridbase, $index;
  system("gridslice k $index <$gridname >$basename.asc 2>/dev/null");
  push(@files, "$basename.asc");

  if ($tilecount % 4 == 0) {
    print STDERR '.';
    if ($tilecount > 0) {
      push(@gnuscript, "unset multiplot");
    }
    my $name = sprintf "%s_sliced_%02d.ps", $gridbase, $tilecount/4;
    push(@gnuscript, "set out \"$name\"");
    push(@gnuscript, "set multiplot layout 2,2");
  }

  push(@gnuscript, "plot \"$basename.asc\" matrix with image title \"$index\"");
  ++$tilecount;
}

my $gnuplot = new FileHandle "|gnuplot 2>/dev/null";
defined($gnuplot) || die "$0: Error- cannot open pipe to gnuplot";
print $gnuplot join("\n", @gnuscript);
undef($gnuplot);

unlink @files;

print STDERR "\nDone.\n";

sub getGridSize {
  my $fn = shift;
  my $fh = new FileHandle "gridstat 5 5 <$fn|";
  defined($fh) || die "$0: Error - cannot open pipe from gridstat command";

  $_ = <$fh>;
  chomp;
  /\((\d+),(\d+),(\d+)\)/;
  return( ($1, $2, $3) );
}


sub parseRanges {
  my $opt = shift;

  my @indices;
  my @ary = split(/,/, $opt);
  foreach (@ary) {
    my($a, $b, $c);
    if (/(\d+):(\d+):(\d+)/) {
      $a = $1;
      $b = $3;
      $c = $2;
    } elsif (/(\d+):(\d+)/) {
      $a = $1;
      $b = $2;
      $c = ($a > $b) ? 1 : -1;
    } else {
      $a = $_;
      $b = $_;
      $c = 1;
    }

    my $i;
    if ($a <= $b) {
      for ($i = $a; $i <= $b; $i += $c) {
	push(@indices, $i);
      }
    } else {
      for ($i = $b; $i >= $a; $i -= $c) {
	push(@indices, $i);
      }
    }
  }

  return(@indices);
}



sub show_help {
print <<EOF;
Usage: grid_slicer.pl [options] gridname
       --stride=i     Sets the stride in k to section the grid
       --range=s      A comma-separated list of Octave-style ranges selecting
                      which k-sections to extract
       --output=s     Sets the output prefix, otherwise the gridname (minus
                      any extensions) is used
       --help         Displays this message
EOF

exit;
}
