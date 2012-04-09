#!/usr/bin/perl
#
# blurrogram - Concatenate structures together to create a composite figure showing motion
#
#
#  This file is part of LOOS.
#
#  LOOS (Lightweight Object-Oriented Structure library)
#  Copyright (c) 2008, 2009 Tod D. Romo
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


=head1 NAME

blurrogram - concatenates structures together from a trajectory to create a composite figure

=head1 SYNOPSIS

blurrogram.pl [options] modelname trajectoryname >script.pml


=head1 DESCRIPTION

blurrogram extracts structures from a trajectory and writes them out
as a set of PDB files.  It also generates a Pymol script that can be
run from within Pymol to load these structures and display them.  This
is useful for showing the range of motion of a structure in a single
figure. 

The Pymol script is composed of three parts: a preamble, the body, and
a postscript.  Each part can be defined separately.  The body of the
script is special in that certain variables are injected into it.
These are "pdbname" for the name of the PDB file, "name" for the
basename (i.e. pdbname without the ".pdb" suffix), iframe" for the
integer frame number from the trajectory and "cox" for the current
output index.

If any component of the script is not specified, a default one will be
inserted.  The default preamble is to turn orthscopic view on and set
the helix length to something small (if the "cartoon" style is
used). The default script body loads each structure and displays it as
a cartoon.  It is then colored using a spectrum.  There is no default
postscript.

=head2 OPTIONS

=over

=item --frames

Sets the number of structures generated

=item --script

Reads the script body from the specified filename or specifies commands to use

=item --preamble

Reads the script preamble from the specified filename or specifies commands to use

=item --postscript

Reads the script postscript from the specified filename or specifies commands to use

=item --filespec

This is a printf-formatted string that gives the base frame filename

=item --selection

This is a LOOS-selection that is passed to the frame2pdb program to
generate the PDB frames

=item --range

A Matlab/Octave-style range of trajectory frames to operate over

=item --style

This sets the display style for the default Pymol script

=item -color

This sets the color for each structure in the default Pymol script


=back

=head1 EXAMPLES

 blurrogram.pl --select='segid =~ "BAR."' b2ar.pdb b2ar.dcd >script.pml
 This extracts 10 frames from b2ar.dcd only using the protein (via the --selection)

 blurrogram.pl --style='ribbon' --color='red, $name', --select='segid =~ "BAR."' b2ar.pdb b2ar.dcd >script.pml
 This is as above, but uses a ribbon style with each structure colored red.

 blurrogram.pl --nframes=20 --script='myscript.pml' b2ar.pdb b2ar.dcd >script.pml
 This extracts 20 frames using "myscript.pml" as the body script rather than the default.
 With the following script, it will display each structure as a red ribbon:
   load $pdbname
   hide lines, $name
   show ribbon, $name
   spectrum count, rainbow, $name

 blurrogram.pl --range=0:100:1000 b2ar.pdb b2ar.dcd >script.pml
 This extracs 11 frames (0, 100, 200, ... 1000)

=cut

use FileHandle;
use Getopt::Long;

### Configurable things...
$trajinfo = 'trajinfo';
$frame2pdb = 'frame2pdb';
$color = 'spectrum count, rainbow, $name';
$style = 'ribbon';


### Process command-line options...

$header = join(' ', @ARGV);

$nframes = 10;
@script = ();
@preamble = ();
@postscript = ();
$filespec = 'frame_%03d';
$selection = 'all';
$helpflag = 0;
$fixed_start = -1;
$fixed_end = -1;
$fixed_step = -1;

$result = GetOptions("frames=i" => \$nframes,
		     "script=s" => sub { @script = &readFromFile($_[1]); },
		     "filespec=s" => \$filespec,
		     "selection=s" => \$selection,
		     "preamble=s" => sub { @preamble = &readFromFile($_[1]); },
		     "postscript=s" => sub { @postscript = &readFromFile($_[1]); },
		     "range=s" => sub {
		       my @opts = @_;
		       my @ary = split(/:/, $opts[1]);
		       if ($#ary == 1) {
			 $fixed_start = $ary[0];
			 $fixed_end = $ary[1];
			 $fixed_step = 1;
		       } elsif ($#ary == 2) {
			 $fixed_start = $ary[0];
			 $fixed_end = $ary[2];
			 $fixed_step = $ary[1];
		       } else {
			 die "$0: Error- cannot parse range '$opts[1]'";
		       }
		       $nframes = ($fixed_end - $fixed_start) / $fixed_step + 1;
		     },
		     "color=s" => \$color,
		     "style=s" => \$style,
		     "help" => \$helpflag);

if ($#ARGV != 1 || $helpflag) {
  &showHelp;
  exit;
}
$model_name = shift;
$traj_name = shift;

if ($#preamble < 0) {
  @preamble = ("set orthoscopic, 1");
  if ($style eq 'cartoon') {
    push (@preamble, "set cartoon_oval_length,0.3");
  }
}

if ($#script < 0) {
  push(@script, "load \$pdbname");
  push(@script, "hide lines, \$name");
  push(@script, "show $style, \$name");
  push(@script, "$color");
}

### Default postscript would go here...

if ($#postscript < 0) {
}


############################################
print "# $0 $header\n";
&generateScript(\@preamble);

# Get # of frames...
$nsteps = &getTrajInfo($model_name, $traj_name);
print STDERR "Found $nsteps frames in $traj_name\n";
$stepsize = ($nsteps-1) / ($nframes-1);
if ($fixed_start < 0) {
  $fixed_start = 0;
  $fixed_step = $stepsize;
  $fixed_end = $nsteps;
}

for ($cox=0; $cox<$nframes; ++$cox) {
  $frame = $fixed_start + $cox * $fixed_step;
  $iframe = int($frame);
  $name = sprintf $filespec, $cox;
  $pdbname = $name . '.pdb';
  print STDERR "Processing blurrogram frame #$cox with trajectory frame $iframe...\n";
  system("$frame2pdb --selection '$selection' $model_name $traj_name $iframe >$pdbname");
  &generateScript(\@script);
  print "\n";
}

&generateScript(\@postscript);



sub showHelp {
  print STDERR <<EOF; 
Usage: blurrogram.pl [options] model trajectory
Options:
   --help          This message
   --frames=i      # of frames to generate
   --script=s      Name of Pymol script to insert, otherwise use default script
   --preamble=s    Filename of the preamble for the generated script
   --postscript=s  Filename of the postscript for the generated script
   --filespec=s    Printf-formatted string specifying frame filenames  ($filespec)
   --selection=s   Selection string to pass to frame2pdb ($selection)
   --range=s       MATLAB/Octave style range of trajectory frames to use
   --style=s       Style passed to "show" command for default script  ($style)
   --color=s       Coloring line used in the default script           ($color)
EOF
}


# Call trajinfo to get the number of frames in the trajectory
sub getTrajInfo {
  my $mname = shift;
  my $tname = shift;

  my $fh = new FileHandle "$trajinfo --brief 1 $mname $tname|";
  my $output = <$fh>;
  chomp($output);
  my @ary = split(/\s+/, $output);
  return($ary[1]);
}


# Reads the contents of a file, returning each line in an array
sub readFromFile {
  my $fn = shift;

# Hack to see if what's passed is a file.  If not, split and return it
# as the array...this lets you put Pymol commands on the command-line
  my @lines;
  if (-e $fn) {
    my $fh = new FileHandle $fn;
    defined($fh) || die "$0: Error- cannot open $fn for reading.";
    @lines = <$fh>;
  } else {
    @lines = split(/[\n;]+/, $fn)
  }

  return(@lines);
}


# Generate the script by eval'ing each line so variables can be substituted
sub generateScript {
  my $rary = shift;

  foreach (@$rary) {
    my $line;
    chomp;
    eval("\$line = \"$_\";");
    print "$line\n";
  }
}
