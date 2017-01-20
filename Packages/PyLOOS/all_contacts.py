#!/usr/bin/python

import loos
import loos.pyloos
import sys
import numpy
from os.path import basename, splitext
import argparse

#if len(sys.argv) < 5 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
def fullhelp():
    print """
        Description: all_contacts.py: tracks contacts between pairs of side
                     chains in a protein or nucleic acid.  Returns a
                     matlab-format matrix containing the fraction of frames
                     containing each residue-residue contact, as well as a
                     matrix for each trajectory individually.

                     The threshold for deciding if 2 residues are in contact is
                     whether any atoms from one residue are within 4.0 ang of
                     any atom from the other.  NOTE: PERIODICITY IS NOT
                     RESPECTED, since we're assuming they're both in the same
                     molecule.
          """
class FullHelp(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        kwargs['nargs']=0
        super(FullHelp, self).__init__(option_strings, dest, **kwargs)
    def __call__(self, parser, namespace, values, option_string = None):
        fullhelp()
        parser.print_help()
        setattr(namespace, self.dest, True)
        parser.exit()

parser = argparse.ArgumentParser(description="Track residue-residue contacts")
# Positional arguments
parser.add_argument('system_file', help="System file, e.g. psf, pdb")
parser.add_argument('selection_string', help="Selection string for analysis")
parser.add_argument('out_file', help="Matrix of total contacts")
parser.add_argument('--skip_hydrogen',
                    help="Don't use hydrogens for contact calculation",
                    action='store_true')
parser.add_argument('--skip_backbone',
                    help="Ignore backbone atoms",
                    action='store_true')
parser.add_argument('--fullhelp',
                    help="Print detailed description of all options",
                    action=FullHelp)
parser.add_argument('--cutoff',
                    help="Cutoff distance",
                    default=4.0,
                    type=float)
parser.add_argument('--threshold',
                    help="Number of atom-atom contacts to say residues touch",
                    default=1,
                    type=int)


# This must be the last argument on the command line
parser.add_argument('--traj_files', help="One of more trajectory files",
                    nargs= argparse.REMAINDER,
                    required=True)

args = parser.parse_args()
header = " ".join(sys.argv) + "\n"
print "#", header


system = loos.createSystem(args.system_file)
all_trajs = []
out_names = []
num_trajs = len(args.traj_files)
for t in traj_files:
    traj = loos.pyloos.Trajectory(t, system)
    all_trajs.append(traj)
    t_base = basename(t)
    core, ext = splitext(t_base)
    out_names.append(core + ".dat")

if args.skip_hydrogen:
    no_hydrogens = loos.selectAtoms(system, "!hydrogen")
    target = loos.selectAtoms(no_hydrogens, selection)

target = loos.selectAtoms(system, selection)


residues = target.splitByResidue()

if args.skip_backbone:
    residues = list(map(lambda r: loos.selectAtoms(r, "!backbone"), residues))

frac_contacts = numpy.zeros([len(residues), len(residues), num_trajs],
                            numpy.float)


for traj_id in range(num_trajs):
    traj = all_trajs[traj_id]
    for frame in traj:
        for i in range(len(residues)):
            for j in range(i+1, len(residues)):
                if residues[i].contactWith(args.cutoff, residues[j], args.threshold):
                    frac_contacts[i, j, traj_id] += 1.0
                    frac_contacts[j, i, traj_id] += 1.0
    frac_contacts[:, :, traj_id] /= len(traj)
    numpy.savetxt(out_names[traj_id], frac_contacts[:, :, traj_id],
                  header=header)

average = numpy.add.reduce(frac_contacts, axis=2)
average /= len(traj_files)


numpy.savetxt(out_file, average, header=header)
