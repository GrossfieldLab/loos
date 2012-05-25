#!/usr/bin/perl -w
#
# 
#  This file is part of LOOS.
#
#  LOOS (Lightweight Object-Oriented Structure library)
#  Copyright (c) 2012 Tod D. Romo
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

hoccupancies - determine occupancies from a hydrogen bond state matrix

=head1 SYNOPSIS

hoccupancies.pl [--frameno] matrix.asc >output.asc


=head1 DESCRIPTION

Computes the occupancies (down the column) for hydrogen bonds.  This
tool can be used with the output from hmatrix or hcontacts.  The
latter includes the frame number as the first column in the matrix.
This can be ignored by using the --frameno option.

=cut


use strict;
use FileHandle;
use Getopt::Long;

my $frameno = 0;

my $hdr = $0 . ' ' . join(' ', @ARGV);

my $ok = GetOptions(
		    'frameno!' => \$frameno,
		    'help' => sub { &help; exit; }
		   );

# Process args
if (!$ok) {
  &help;
  die;
}

# Capture existing metadata
my $meta;
my $rows;
my $cols;

while (<>) {
  if (/^# (\d+) (\d+) \(\d\)/) {
    $rows = $1;
    $cols = $2;
    last;
  }

  if (/^#/) {
    $meta .= $_;
  }
}

my $firstcol = ($frameno) ? 1 : 0;

my @counts = (0) x $cols;

while (<>) {
  chomp;
  my @ary = split;
  $#ary+1 == $cols || die "Error- mismatch in number of columns at line $.";
  for (my $i=$firstcol; $i<$cols; ++$i) {
    $counts[$i] += $ary[$i];
  }
}


print "# $hdr\n";
print "$meta";
print "# ", $frameno ? $cols-1 : $cols, " 2 (0)\n";
for (my $i=$firstcol; $i<$cols; ++$i) {
  print $i, "\t", $counts[$i] / $rows, "\n";
}


sub help {
  print "Usage- hoccupancies [--frameno] matrix.asc\n";
}
