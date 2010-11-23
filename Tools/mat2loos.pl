#!/usr/bin/perl -w
#
#  This file is part of LOOS.
#
#  LOOS (Lightweight Object-Oriented Structure library)
#  Copyright (c) 2010 Tod D. Romo, Grossfield Lab
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
use File::Temp qw/tempfile tempdir/;


my $replace = 1;
my $stdout = 0;

my $tmpdir = '.';
my $template = 'matrix_XXXXXX';
my $tmpsuffix = '.asc';


my %opts = (
	    'stdout' => sub { $replace=0; $stdout=1; },
	    'replace' => sub { $replace=1; $stdout=0; },
	    'help' => sub {&help;}
	   );

$hdr = "# $0 " . join(' ', @ARGV);

GetOptions(%opts) || &help;
$#ARGV == 0 || &help;


my $filename = shift;
my($m, $n) = &findSize($filename);

my $fin = new FileHandle $filename; defined($fin) || die "Error- cannot open $filename";
my $fout = $replace ? File::Temp->new(TEMPLATE => $template, DIR => $tmpdir, SUFFIX => $tmpsuffix, UNLINK => 0) : STDOUT;
defined($fout) || die "Error- cannot open output file";

print $fout $hdr, "\n";
print $fout "# $m $n (0)\n";
while (<$fin>) {
  print $fout $_;
}

if ($replace) {
  unlink $filename;
  rename $fout->filename, $filename;
}



sub findSize {
  my $fn = shift;
  my $fh = new FileHandle $fn; defined($fh) || die "Error- cannot open $fn";

  my $m = 0;
  my $n = 0;

  while (<$fh>) {
    next if /^#/;
    if ($n == 0) {
      my @ary = split;
      $n = $#ary+1;
    }
    ++$m;
  }

  return($m, $n);
}



sub help {
  print <<EOF;
Usage- mat2loos.pl [--replace|--stdout] matrix.asc

   --replace   Replaces the given file with the new LOOS-formatted matrix
   --stdout    Sends the new LOOS-formatted matrix to stdout

NOTES:
 o Does not work on triangular or sparse matrices
EOF

  exit 0;
}
