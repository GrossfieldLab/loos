#!/usr/bin/perl -w
#
# (c) 2011 Tod D. Romo, Grossfield Lab, URMC
#
# An attempt to verify that a LOOS install is "correct", or at least
# not missing anything...
#
# Usage- check_loos_install.pl  loos_source_dir    loos_install_dir



use strict;
use FileHandle;
use File::Basename;

use Getopt::Long;

my $full = 1;
my $prep = 0;
my @excludes = ();

my $ok = GetOptions(
		    'full!' => \$full,
		    'prep!' => \$prep,
		    'exclude=s' => \@excludes
		   );

$ok || die "Error- usage check_loos_install.pl [--[no]full] [--[no]prep] source_dir install_dir";

my $srcdir = shift;
my $trgdir = shift;


my $rsrcfiles;
my $rtrgfiles;
my $n;
($rsrcfiles,$n) = &findFiles($srcdir);
print "Read in $n files from $srcdir\n";

($rtrgfiles,$n) = &findFiles($trgdir);
print "Read in $n files from $trgdir\n";

my %missing;
my $nmissing = 0;
foreach my $filename (keys %$rsrcfiles) {
  next if ($rsrcfiles->{$filename}->{ADSUM});
  if ($rtrgfiles->{$filename}) {
    $rsrcfiles->{$filename}->{ADSUM} = 1;
  } else {
    ++$nmissing;
    my $path = $rsrcfiles->{$filename}->{PATH};
    if (!exists($missing{$path})) {
      $missing{$path} = [ $filename ];
    } else {
      push(@{$missing{$path}}, $filename);
    }
  }
}

if ($nmissing > 0) {
  print "***There were $nmissing missing files***\n";
  foreach my $path (sort keys %missing) {
    my @filelist = @{$missing{$path}};
    print "$path:\n";
    foreach my $filename (@filelist) {
      print "\t$filename\n";
    }
  }
  exit -1;
}

exit 0;



sub findFiles {
  my $dir = shift;

  my $n = 0;
  my $find = new FileHandle "find $dir \\( -name '*.hpp' -o -name '*.h' -o \\( -perm /111 -a ! -type d \\) \\) -print|";
  defined($find) || die "Error- cannot open pipe from find command for $dir";

  my %files;
  while (<$find>) {
    next if /.sconf/;
    if (! $full) {
      next if /Packages/;
    }
    if ($prep) {
      next if /prep_release.sh/;   # Hack!
    }

    my $skipit = 0;
    foreach my $exclude (@excludes) {
      if (/$exclude/) {
	$skipit = 1;
	last;
      }
    }

    next if $skipit;

    chomp;
    my($name, $path, undef) = fileparse($_);
    $files{$name} = { NAME => $name, PATH => $path };
    ++$n;
  }

  return(\%files, $n);
}
