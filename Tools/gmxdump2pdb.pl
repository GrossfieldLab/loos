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



#############################################################################
#
# Translate a subset of a TPR file (via gmxdump) into a PDB file
# complete with CONECT records and a fake-PSF
#
# Usage- gmxdump -s foo.tpr | gmxdump2pdb.pl prefix
#
#  This generates prefix.pdb and prefix.psf
#
# NOTES:
#
#  o This program was written by observing the output of gmxdump for
#    our test-cases.  No guarantee is made that it will work for any
#    specific Gromacs/MARTINI system.
#
#  o The PSF is a fake-PSF...it only contains atom types and
#    bonds, but it will work with VMD and LOOS
#
#  o The PDB will contain CONECT records that also describe
#    the bonds in the system.  These can be read into pymol
#    by first typing "set connect_mode, 1" in pymol, and then
#    loading the PDB.  These are also usable in LOOS (i.e.
#    you can use the PDB file to define both coordinates and
#    connectivity rather than using a PSF.  For large systems,
#    however, the PSF will be read faster by LOOS.)
#
#  o Residue id's (internal to molecules) are taken from the bracketed
#    number in the topology, not from the "nr=\d" field...
#
#  o Hydrogens may appear as constraints.  Add them in with
#    either --constraints (to get all constraints) or
#    --hydrogens to get only the hydrogens.  Note that water
#    is treated differently and bonds will not be made with just --hydrogens
#
#  o Water h-bonds can be inferred using the --water option (requires
#    water to consist of OW, HW1, and HW2 atoms in that order
#
#############################################################################


use FileHandle;
use Carp;
use Getopt::Long;
use strict;
use POSIX qw(ceil);

my $debug = 0;    # Enable copious output for debugging...

my $resids_local = 0;
my $coord_scale = 10.0;   # Convert NM into Angstroms for PDBs
my $use_constraints = 0;
my $use_vsites = 0;
my $hydrogens_only = 0;
my $infer_water = 0;
my $force_mass = 0;
my $tabbonds = 1;
my $conbonds = 1;
my $namd = 1;


my %seen_segids;     # Mapping of gromacs molecule segment names to SEGIDs
my %used_segids;     # Also keep track of user-specified segid mappings to make sure
                     # all were used


my $ok = GetOptions('local!' => \$resids_local,
                    'constraints!' => \$use_constraints,
                    'vsites!' => \$use_vsites,
                    'hydrogens' => sub { $use_constraints = $hydrogens_only = 1; },
                    'water!' => \$infer_water,
                    'mass=f' => \$force_mass,
                    'tabbonds!' => \$tabbonds,
                    'conbonds!' => \$conbonds,
                    'segmap=s' => \%seen_segids,
                    'namd!' => \$namd,
                    'help' => sub { &showHelp; });

$ok || &showHelp;


# Mark user-supplied segid mappings
foreach (keys(%seen_segids)) {
  $used_segids{$_} = 0;
  my $val = $seen_segids{$_};
  if (length($val) > 4) {
    die "Error- you cannot use a segid that is longer than 4 characters.";
  }
}

my $bondprefix_regex = &build_bond_prefix_regex($tabbonds, $conbonds);

my $output_prefix = shift(@ARGV);

my $pdb = new FileHandle "$output_prefix.pdb", 'w'; defined($pdb) || die "Error- cannot open $output_prefix.pdb for writing";
my $psf = new FileHandle "$output_prefix.psf", 'w'; defined($psf) || die "Error- cannot open $output_prefix.psf for writing";

# Processing the gmxdump output is done in 3 stages.  First, the
# overall structure is parsed.  Then the detailed information about
# each molecule is filled in.  Finally, the coordintes are read in.

my $rtopo = &processTopology;
print STDERR "Found ", $#$rtopo+1, " segments in topology\n";

my $rmolecules = &processMolecules;
print STDERR "Found ", $#$rmolecules+1, " segments defined\n";

my $rbox = &processBox;

my $rcoords = &readCoords;
print STDERR "Found ", $#$rcoords+1, " coordinates.\n";

# Here, we build up an array of ATOMs roughly with PDB fields, along
# with the appropriate connectivity information.

my($rstruct, $rconn) = &buildStructure($rtopo, $rmolecules, $rcoords);
&inferWater($rstruct, $rconn) if ($infer_water);

# Write out the PDB.  Uses hybrid-36 format if resid or atomid overflows...
print $pdb "REMARK    MADE BY GMXDUMP2PDB.PL\n";
printf $pdb "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%11s%4d\n",
  10.*$rbox->[0]->[0], 10.*$rbox->[1]->[1], 10.*$rbox->[2]->[2],
  90.0, 90.0, 90.0, 'P1', 1;

my $natoms = 0;

foreach my $atom (@$rstruct) {
  # Try to catch missing information...
  defined($$atom{ATOMID}) || die;
  defined($$atom{ATOMNAME}) || die;
  defined($$atom{X}) || die;
  defined($$atom{RESNAME}) || die;
  defined($$atom{RESID}) || die;
  defined($$atom{SEGID}) || die;

  my $line = sprintf("%-6s%5s %4s %4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n",
		     'ATOM',
		     $$atom{ATOMID} >= 100000 ? &toHybrid36($$atom{ATOMID}, 5) : $$atom{ATOMID},
		     $$atom{ATOMNAME},
		     $$atom{RESNAME},
		     ' ',
		     $$atom{RESID} >= 10000 ? &toHybrid36($$atom{RESID}, 4) : $$atom{RESID},
		     $$atom{X},
		     $$atom{Y},
		     $$atom{Z},
		     1.0,
		     1.0,
		     &sanitizeSegid($$atom{SEGID}));
  print $pdb $line;
  ++$natoms;

}

my $expected_natoms = $#$rstruct + 1;
$natoms == $expected_natoms || warn "Warning- wrote $natoms, but expected $expected_natoms\n";

my @bond_list = sort { $a <=> $b } keys %$rconn;
my $total_bonds = 0;

foreach my $atom (@bond_list) {

  # Bound atoms may appear multiple times, especially if adding in
  # constraints.  Filter out the duplicates...
  my $rbound = &uniqueElements($$rconn{$atom});
  my $n = 0;

  foreach my $bound (@$rbound) {
    ++$total_bonds;

    printf $pdb "CONECT%5s", ($atom >= 100000 ? &toHybrid36($atom, 5) : $atom) if ($n == 0);
    printf $pdb '%5s', $bound >= 100000 ? &toHybrid36($bound, 5) : $bound;
    if (++$n > 5) {
      print $pdb "\n";
      $n = 0;
    }
  }
  print $pdb "\n" if ($n != 0);
}


########
######## Build a fake PSF
########

# Note: This is the charmm extended PSF format.  There is an issue with VMD
#       where some fields are truncated by 1 char (i.e. resname is 7 chars rather than
#       8).  We get around this by left-justifying all fields.  But the field is still
#       smaller than expected...

print $psf "PSF EXT", ($namd ? ' NAMD' : ''), "\n\n";
print $psf "       2 !NTITLE\nREMARKS THIS IS NOT A REAL PSF, USE CAREFULLY\nMADE BY GMXDUMP2PDB.PL\n\n";
printf $psf "%8d !NATOM\n", $natoms;
foreach my $atom (@$rstruct) {
  defined($$atom{ATOMID}) || die "Undefined atomid";
  defined($$atom{SEGID}) || die "Undefined segid";
  defined($$atom{RESID}) || die "Undefined resid";
  defined($$atom{RESNAME}) || die "Undefined resname";
  defined($$atom{ATOMNAME}) || die "Undefined atomname";
  defined($$atom{CHARGE}) || die "Undefined charge for $$atom{ATOMID}";
  defined($$atom{MASS}) || die "Undefined mass";

  my $mass = $force_mass>0.0 ? $force_mass : $$atom{MASS};

  printf $psf "%-10s %-8s %-8s %-8s %-8s %-4s %-14.6g%-14.6g%8d\n",
    $$atom{ATOMID},
      &sanitizeSegid($$atom{SEGID}),
	$$atom{RESID},
	  $$atom{RESNAME},
	    $$atom{ATOMNAME},
	      $$atom{ATOMTYPE},
		$$atom{CHARGE},
		  $mass,
		    0;

}

# Now right out the bond pair list...
printf $psf "\n%8d !NBOND: bonds\n", $total_bonds;
my $count = 0;
foreach my $atom (@bond_list) {
  my $rbound = &uniqueElements($$rconn{$atom});
  foreach (@$rbound) {
    printf $psf "%10s%10s", $atom, $_;
    if (++$count >= 4) {
      $count = 0;
      print $psf "\n";
    }
  }
}
print $psf "\n";


&warnIfUnusedSegmentMaps;

# Returns a reference to an array containing references to hashes that
# represent each molecule contained in the TPR.  It is assumed that
# "ffparams:" in the output denotes the likely end of this block of
# data...
#
# Tags:
#   NAME -> name of the molecule (i.e. segid)
#   TYPE -> moltype
#   NMOLS -> Number of molecules present (for this molecule)
#   NATOMS -> Number of atoms in each molecule

sub processTopology {
  my @molecules;     # Array of anon-hashes containing molecule info
  my $state = 0;     # Initial state for our state machine...


  my $curmol = undef;   # The molecule we're currently parsing

  while (<>) {
    if ($state == 0) {

      if (/^topology:/) {
	$state = 1
      }

    } elsif ($state == 1) {

      if (/moltype\s+= (\d+) \"(.+)\"/) {
	if (defined($curmol)) {
	  push(@molecules, $curmol);
	  $curmol = {};
	}

	$curmol->{NAME} = $2;
	$curmol->{TYPE} = $1;
	print STDERR "Found segment '$2' [type=$1] (", $#molecules+1, ")\n";
      } elsif (/molecules\s+= (\d+)/) {
	$curmol->{NMOLS} = $1;
	print STDERR "\tMolecules: $1\n";
      } elsif (/atoms_mol\s+= (\d+)/) {
	$curmol->{NATOMS} = $1;
	print STDERR "\tAtoms: $1\n";
      } elsif (/ffparams:/) {
	push(@molecules, $curmol);
	last;
      }

    }

  }

  return(\@molecules);
}


# Process the rest of the gmxdump output, adding the appropriate atoms, coords, and bonds.
#
# We assume "grp[" denotes the end of the metadata block (at least as
# much of it as we care about).
#
# Added tags:
#   ATOMS  ->  anon-array of atom index to name
#   CHARGES -> anon-array of atom index to charge
#   TYPES -> anon-array of atom index to type
#   MASSES -> anon-array of atom index to mass
#   RESIDUES -> anon-hash of residue index to name
#   ATOM_TO_RESIDUE -> anon-hash of indices in residues for corresponding atom
#   BONDS -> anon-array of pairs of indices to atoms that are bound
#   CONSTRAINTS -> anon-array to constraints (to be used as bonds)
#   VSITES -> anon-array of virtual site constraints to be used as bonds
#             (note: only vsite type 2 is supported)


sub processMolecules {
  my $rmols = [];       # Molecules we'll be filling in info for...
  # the details for...
  my $state = 0;        # Initial state for our state machine

  my $molidx = -1;   # Index into rmols for what we're processing
  my $ratoms = undef;   # Anon-array of atom names for each molecule
  # as we're processing it
  my $rbonds = undef;   # Same, but for connectivity
  my $rcons = undef;    # Constraints
  my $rvsites = undef;



  while (<>) {
    print "[$state : $molidx] $_" if ($debug);

    if ($state == 0) {

      if (/moltype \((\d+)\)/) {
        $molidx = $1;
        $state = 1;
      } elsif (/grp/) {
        last;
      }

    } elsif ($state == 1) {

      if (/name="(.+)"/) {
      }

      # Found an atom info section, so instantiate a new array to hold
      # the atom names and link it in...  Since there are two atom
      # blocks in the output, and only the second one is the one we're
      # interested in, delay actually processing it until the second
      # (via state 3)
      if (/atom \(\d+\)/) {
        $state = 2;
        $$rmols[$molidx]->{ATOMS} = [];
        $$rmols[$molidx]->{CHARGES} = [];
        $$rmols[$molidx]->{TYPES} = [];
        $$rmols[$molidx]->{MASSES} = [];

        $$rmols[$molidx]->{ATOM_TO_RESIDUE} = {};
        $$rmols[$molidx]->{RESIDUES} = {};
      }

    } elsif ($state == 2) {

      # We've started the 2nd block of atoms, so process subsequent
      # atom lines
      if (/atom \(\d+\)/) {
        $state = 3;
      } elsif (/atom\[\s*(\d+)\].*res(ind|nr)=\s*(\d+)/) {
        $$rmols[$molidx]->{ATOM_TO_RESIDUE}->{$1} = $3;
        if (/ q=([ 0-9.eE+-]+),/) {
          push(@{$$rmols[$molidx]->{CHARGES}}, $1 * 1.0);   # force to be float
        }
        if (/ m=([ 0-9.eE+-]+),/) {
          push(@{$$rmols[$molidx]->{MASSES}}, $1 * 1.0);    # force to be float
        }
      }

    } elsif ($state == 3) {

      if (/moltype \((\d+)\)/) {
        $molidx = $1;
        $state = 1;

      } elsif (/atom\[(\d+)\]=\{name="(.+)"\}/) {
        $$rmols[$molidx]->{ATOMS}->[$1] = $2;

      } elsif (/residue\[(\d+)\]=\{name="(.+)"/) {
        # NOTE: This ignores the "nr=\d+" field in lieu of the
        # bracketed index.  I'm not sure if this is the correct thing
        # to do...
        $$rmols[$molidx]->{RESIDUES}->{$1} = $2;

      } elsif (/type\[(\d+)\]=\{name="([^"]+)",/) {
        push(@{$$rmols[$molidx]->{TYPES}}, $2);

      } elsif (/Bonds?:/) {
        $state = 4;
        $rbonds = [];
        $$rmols[$molidx]->{BONDS} = $rbonds;
      } elsif (/moltype \((\d+)\)/)  {
        $molidx = $1;
        $state = 1;
      } elsif ($use_constraints && /Constraint/) {
        $rcons = [];
        $$rmols[$molidx]->{CONSTRAINTS} = $rcons;
        $state = 5;
      } elsif (/groupnr\[/) {
        last;
      }

    } elsif ($state == 4) {
      if (/\(($bondprefix_regex)?BONDS\)\s+(\d+)\s+(\d+)/) {
        my $rpair = [$2, $3];
        push(@$rbonds, $rpair);
      } elsif ($use_constraints && /Constraint/) {
        $rcons = [];
        $$rmols[$molidx]->{CONSTRAINTS} = $rcons;
        $state = 5;
      } elsif ($use_vsites && /Virtual site (\d):/) {
        if ( $1 != 2 ) {
          warn "Warning- found unsupported virtual site type $1 in topology.  Will skip...";
        } else {
          $rvsites = [];
          $$rmols[$molidx]->{VSITES} = $rvsites;

          $state = 7;
        }
      } elsif (/moltype \((\d+)\)/)  {
        $molidx = $1;
        $state = 1;
      } elsif (/groupnr\[/) {
        last;
      }

    } elsif ($state == 5) {
      if (/\(CONSTR\) (\d+) (\d+)/) {
        my $rpair = [$1, $2];
        push(@$rcons, $rpair);
      } elsif ($use_vsites && /Virtual site (\d)(out)?:/) {
        if ( $1 != 2 && $1 != 3 ) {
          warn "Warning- found unsupported virtual site type $1 in topology.  Will skip...";
        } else {
          $rvsites = [];
          $$rmols[$molidx]->{VSITES} = $rvsites;

          $state = 7;
        }
      } elsif (/moltype \((\d+)\)/) {
        $molidx = $1;
        $state = 1;
      } elsif (/groupnr\[/) {
        last;
      } elsif (/Constr. No Conn./) {
        #	$state = 6;
      }

      ### Disabled for now...
    } elsif ($state == 6) {
      if (/\(CONSTRNC\) (\d+) (\d+)/) {
        my $rpair = [$1, $2];
        push(@$rcons, $rpair);
      } elsif ($use_vsites && /Virtual site (\d):/) {
        if ($1 != 2) {
          warn "Warning- found unsupported virtual site type $1 in topology.  Will skip...";
        } else {
          $rvsites = [];
          $$rmols[$molidx]->{VSITES} = $rvsites;

          $state = 7;
        }
      } elsif (/moltype \((\d+)\)/) {
        $molidx = $1;
        $state = 1;
      } elsif (/groupnr\[/) {
        last;
      }

    } elsif ($state == 7) {
      if (/\(VSITE(\d)(OUT)?\) ([0-9 ]+)$/) {
        ($1 == 2 || $1 == 3) || die "Error- cannot process vsite type $1";
        my @indices = split(/\s+/, $3);
        if ($#indices == 2) {
          push(@$rvsites, [$indices[0], $indices[1]]);
          push(@$rvsites, [$indices[1], $indices[2]]);
        } elsif ($#indices == 3) {
          push(@$rvsites, [$indices[0], $indices[1]]);
        }
      } elsif (/moltype \((\d+)\)/) {
        $molidx = $1;
        $state = 1;
      } elsif (/groupnr\[/) {
        last;
      }

    } else {
      croak "Error- unknown state #$state";
    }
  }

  return($rmols);
}


sub processBox {
  my @box;

  while (<>) {
    if (/^box\s+\((\d+)x(\d+)\):/) {
      if ($1 ne '3' || $2 ne '3') {
	warn 'Warning- only box size handled is 3x3.  Current box will be ignored.';
      } else {
	for (my $i = 0; $i < 3; ++$i) {
	  $_ = <>;
	  chomp;
	  if (/box\[\s*(\d+)]=\{(.+)\}/) {
	    my $idx = $1;
	    ($idx >= 0 && $idx <= 2) || die "Error- bad index ($idx) for periodic box vector";
	    my @ary = split(/,/, $2);
	    $box[$idx] = \@ary;
	  }
	}
	last;
      }
    }
  }
  return(\@box);
}


## Read in coordinates as 3-tuples
## Returns ref to array of 3-tuples
sub readCoords {
  my @coords;

  while (<>) {
    chomp;
    if (/^\s+x\[\s*\d+\]=\{(.+)\}$/) {
      my $c = $1;
      $c =~ s/\s//g;
      my @coord = split(/,/, $c);
      push(@coords, \@coord);
    }
  }

  return(\@coords);
}



## Builds up the structure based on the topology, molecule, and coord info.
## Returns a ref to an array of atoms (anon-hashes) and a ref to a
## hash of bonds.  Bonds are keyed by atomid and each entry points to
## an anon-array that contains the atomid's of bound atoms.

sub buildStructure {
  my $rtopo = shift;
  my $rmolecules = shift;
  my $rcoords = shift;

  my @atoms;    # Array of anon-hashes containing ATOM record info
  my %bonds;    # Hash of bonds, keyed by atomid, containing
  # anon-arrays of atomids bound to the key.
  my $atomid = 1;  # Global (within the PDB, that is) atomid
  my $resid = 1;   # Global (within the PDB)  This means that the output resids do not
  # match what GROMACS was using

  foreach my $segment (@$rtopo) {
    my $name = $$segment{NAME};
    my $type = $$segment{TYPE};

    my $mol = $$rmolecules[$type];

    my $x = $mol->{ATOMS};
    if ($#$x < 0) {
      croak "$name [$type] has no atoms";
    }

    if (!exists($segment->{NATOMS})) {
	warn "Warning: Missing atom count in topology so inferring number of atoms in segment $name";
	$segment->{NATOMS} = $#{$mol->{ATOMS}} + 1;
    }

    my $residues = $mol->{RESIDUES};
    my $charges = $mol->{CHARGES};

    my $masses = $mol->{MASSES};
    my $atomtypes = $mol->{TYPES};
    my $atom_to_residue = $mol->{ATOM_TO_RESIDUE};
    my $segid = $$segment{NAME};

    for (my $j=0; $j<$segment->{NMOLS}; ++$j) {
      print "Processing segment $segid [$j]\n" if ($debug);
      my @localids;     # Track atomids within a molecule for creating
      # the connectivity info...
      for (my $i=0; $i<$segment->{NATOMS}; ++$i) {
        my %atom;

        $atom{ATOMID} = $atomid;

        $atom{ATOMNAME} = $mol->{ATOMS}->[$i];
        $atom{SEGID} = $segid;
        $atom{CHARGE} = $charges->[$i];
        defined($atom{CHARGE}) || die "Missing charge for atom $i in molecule $j of block $$mol{NAME}";

        $atom{MASS} = $masses->[$i];
        defined($atom{MASS}) || die "Missing mass for atom $i in molecule $j of block $$mol{NAME}";

        $atom{ATOMTYPE} = $atomtypes->[$i];
        push(@localids, $atomid);

        # Figure out the local resid using the molecule's mapping of atoms to residue
        my $residx = $atom_to_residue->{$i};
        defined($residx) || croak "Error- cannot find residue index for atom $i in mol $segid ($j)";
        $atom{RESNAME} = $residues->{$residx};

        # Apply a global resid

        my $local_resid = $resids_local ? $residx+1 : $resid + $residx;
        $atom{RESID} = $local_resid;

        print "\tProcessing $atomid - $atom{ATOMNAME} - $atom{RESNAME} : local_resid=$local_resid, resid=$resid, residx=$residx\n" if ($debug);

        my $c = $rcoords->[$atomid-1];
        $atom{X} = $$c[0] * $coord_scale;
        $atom{Y} = $$c[1] * $coord_scale;
        $atom{Z} = $$c[2] * $coord_scale;

        push(@atoms, \%atom);

        ++$atomid;
      }

      # Need to increment the global resid counter by the number of residues
      # in this molecule.  Cache the size so we don't have to keep computing it...
      if (! exists $mol->{NRESIDUES}) {
        my @ary = keys %$residues;
        $mol->{NRESIDUES} = $#ary+1;
      }
      $resid += $mol->{NRESIDUES};

      # Now add bonds...
      my $rbonds = $mol->{BONDS};
      if ($use_constraints && defined($mol->{CONSTRAINTS})) {
        # Only create bonds to hydrogens...
        if ($hydrogens_only) {
          foreach my $rpair (@{$mol->{CONSTRAINTS}}) {
            my($a, $b) = @$rpair;
            my $id = $localids[$b];
            if ($atoms[$b]->{ATOMNAME} =~ /^H/) {
              push(@$rbonds, $rpair);
            }
          }
        } else {
          push(@$rbonds, @{$mol->{CONSTRAINTS}});
        }
      }

      if ($use_vsites && defined($mol->{VSITES})) {
        push(@$rbonds, @{$mol->{VSITES}});
      }

      for (my $i=0; $i<=$#$rbonds; ++$i) {
        my($a, $b) = @{$$rbonds[$i]};
        $a = $localids[$a];
        $b = $localids[$b];

        print "Adding bond from $a to $b\n" if ($debug);

        if (!exists($bonds{$a})) {
          $bonds{$a} = [ $b ];
        } else {
          push(@{$bonds{$a}}, $b);
        }
      }
    }
  }

  return(\@atoms, \%bonds);
}





## Basic hybrid-36 support (translated from LOOS implementation)

sub toHybrid36 {
  my $d = shift;
  my $n = shift;

  my $n10 = 10**$n;
  my $n36 = 36**($n-1);

  my $cuta = $n10+$n36*26;
  my $negative = 0;

  if ($d < 0 ) {
    $negative = 1;
    $d = -$d;
  }

  ($d < $n10+52*$n36)  || confess "Error- $d is out of range";

  my $coffset = ord('0');
  my $ibase = 10;

  if ($d >= $cuta) {
    $coffset = ord('a') - 10;
    $ibase = 36;
    $d -= $cuta;
    $d += 10*$n36;
  } elsif ($d >= $n10) {
    $coffset = ord('A') - 10;
    $d -= $n10;
    $d += 10*$n36;
    $ibase = 36;
  }

  my @ary;
  while ($d > 0) {
    my $digit = $d % $ibase;
    $digit += ($digit > 9) ? $coffset : ord('0');
    push(@ary, chr($digit));
    $d = int($d/$ibase);
  }

  if ($negative) {
    push(@ary, '-');
  }

  for (my $i=$#ary+1; $i<$n; ++$i) {
    push(@ary, ' ');
  }

  my $result = '';
  for (my $i = $#ary; $i>=0 ; --$i) {
    $result .= $ary[$i];
  }

  return($result);
}



sub uniqueElements {
  my $rv = shift;
  my %elements;

  foreach (@$rv) {
    $elements{$_} = 1;
  }

  my @unique = keys(%elements);
  return(\@unique);
}


sub inferWater {
  my $rstruct = shift;
  my $rconn = shift;

  for (my $i=0; $i<$#$rstruct-2; ++$i) {
    if ($$rstruct[$i]->{ATOMNAME} eq 'OW' &&
	$$rstruct[$i+1]->{ATOMNAME} eq 'HW1' &&
	$$rstruct[$i+2]->{ATOMNAME} eq 'HW2') {

      my $ow = $$rstruct[$i]->{ATOMID};
      my $hw1 = $$rstruct[$i+1]->{ATOMID};
      my $hw2 = $$rstruct[$i+2]->{ATOMID};

      if (!exists($$rconn{$ow})) {
	$$rconn{$ow} = [$hw1, $hw2];
      } else {
	push(@{$$rconn{$ow}}, ($hw1, $hw2));
      }

    }
  }
}



sub sanitizeSegid {
  my $seg = shift;

  $seg =~ y/a-z/A-Z/;

  if (exists($seen_segids{$seg})) {
    $used_segids{$seg} = 1;
    return($seen_segids{$seg});
  }

  my $oseg = $seg;

  # Try squeezing out underscores
  if (length($seg) > 4) {
    $seg =~ s/_//g;
  }

  # Remove non-alphanum
  $seg =~ s/[^a-zA-Z0-9]//;

  # Preserve first and last letters, try to find suitable middle 2 chars...

  if (length($seg)>4) {
    # letter frequency for English text
    my @letfreqs = qw/E I A N S O R T L C U D P M H G Y B F V K W Z X J Q/;

    my @arg = split(//, $seg);


    my $first = shift(@arg);
    my $lastl = pop(@arg);


    while ($#arg > 1) {
      my $l = shift(@letfreqs);
      last if !defined($l);    # Out of letters

      my $flag;
      do {
	$flag = 0;
	for (my $i=0; $i<=$#arg; ++$i) {
	  if ($arg[$i] eq $l) {
	    splice(@arg, $i, 1);
	    $flag = 1;
	    last;
	  }
	}
      } while ($flag && $#arg > 1);
    }

    $seg = $first . join('', @arg) . $lastl;
  }

  # If it's still too long, just take the first 4 chars
  if (length($seg) > 4) {
    $seg = substr($seg, 0, 4);
  }

  if ($seg ne $oseg && exists($seen_segids{$seg})) {
    my $i;
    $seg = substr($seg, 0, 3);
    for ($i = 0; $i <= 9; ++$i) {
      my $s = $seg . $i;
      if (!exists($seen_segids{$s})) {
	$seg = $s;
	last;
      }
    }
    if ($i >= 10) {
      warn "Warning- some long segids could not be mapped to short equivalents...";
      $seg = substr($oseg, -4);
    }
  }

  $seen_segids{$oseg} = $seg;
  $used_segids{$oseg} = 1;

  if ($oseg ne $seg) {
    print STDERR "Warning- Automatically mapping segment '$oseg' to '$seg'\n";
  }
  return($seg);
}



sub warnIfUnusedSegmentMaps {
  foreach (sort(keys(%used_segids))) {
    if (!($used_segids{$_})) {
      print STDERR "Warning- segid mapping for segment $_ was not used.\n";
    }
  }
}

sub build_bond_prefix_regex {
  my($tabbond, $conbond) = @_;
  my @tags;

  if ($tabbond) {
    push(@tags, 'TAB');
  }
  if ($conbond) {
    push(@tags, 'CONN');
  }

  my $reg = join('|', @tags);

  return($reg);
}



sub showHelp {
  print <<'EOF';
Usage- gmxdump -s foo.tpr | gmxdump2pdb.pl [options] output-file-prefix

Options:
   --[no]local       Resids are local to each molecule (default is off)
                     (i.e. reset based on GROMACS' notion of a molecule)
   --[no]constraints Add constraints as bonds (default is off)
   --[no]vsites      Add virtual sites as bonds (default is off) {Only type 2 supported}
   --hydrogens       Only add constaints where the 2nd atom begins with an H
   --water           Infer water connectivity (requires OW, HW1, and HW2 atoms in order)
   --mass=x          Set all atom masses to x
   --[no]tabbonds    Include TABBONDS in bonds list (default is on)
   --[no]conbonds    Include CONNBONDS in bonds list (Gromacs 5+, default is on)
   --[no]namd        Add NAMD tag to PSF (may help with VMD)
   --segmap key=val  Manually specify a mapping of Gromacs molecule names to PDB segids.
                     The key is the gromacs name and "val" is the resulting segid.
                     The key is also case insensitive.  This option may be used multiple
                     times on the command line to map several different names,
                       gmxdump2pdb.pl --segmap FENGYCIN=FENG --segmap NA+=SOD
EOF

  exit 0;
}
