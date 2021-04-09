#!/usr/bin/env python3

"""
Track a base stacking through a trajectory.  Intended for use with
nucleic acids, to identify base stacking.
"""

"""

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2021 Alan Grossfield, Grossfield Lab
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
import sys
import loos
import loos.pyloos
import numpy as np
import argparse
import json


def fullhelp():
    print("""
    all_stacking.py tracks the stacking of a set of cylindrical objects. In
    principle, these could be anything, but the program is designed to work
    with nucleic acid bases, where it will naturally track the average
    stacking score for each base in the selection with all others.

    Usage: all_stacking.py [OPTIONS] system_file selection_string outfile traj1 [traj2 ...]
        -- system_file is a PDB, PSF, prmtop, etc
        -- selection_string identifies the atoms to be analyzed. They will be
           split by residue number to identify the individual sets to be
           evaluated.
        -- outfile will be a numpy-style text-file containing the results. It
           will be a square matrix (num_res x num_res, where num_res is the
           number of residues in the selection)
        -- traj1 is a trajectory file, eg DCD, xtc. Its contents must precisely
           match the system file. Optionally, you can include more than 1
           trajectory file


    Options
        --skip: # of residues to exclude from the front of each trajectory
        --stride: how to step through the trajectory (eg --stride 10 will read
              every 10th frame)
        --fullhelp: print this message
        --base_carbons: apply an additional selection string to the one given
            on the command line, which will exclude the backbone and only
            include carbons atoms. This is for convenience only, since you
            could do this yourself with the command line selection. We exclude
            hydrogens because for RNA the base hydrogens are generally coplanar
            with the carbons, and we exclude the nitrogens because in some
            bases they're out of the plane.
        """)


class FullHelp(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        kwargs['nargs'] = 0
        super(FullHelp, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        fullhelp()
        parser.print_help()
        setattr(namespace, self.dest, True)
        parser.exit()


parser = argparse.ArgumentParser(description="Track stacking")
parser.add_argument('system_file', help="File describing the system")
parser.add_argument('selection_string',
                    help="Selection string describing which residues to use")
parser.add_argument('outfile_name',
                    help="File with the average contact occupancies")
parser.add_argument('traj_files', nargs='+')
parser.add_argument('--skip', type=int, default=0,
                    help="# of frames to skip")
parser.add_argument('--stride', type=int, default=1,
                    help="Ready every nth frame")
parser.add_argument('--base_carbons', action='store_true',
                    help="Just select the carbons in the base")
parser.add_argument('--fullhelp',
                    help="Print detailed description of all options",
                    action=FullHelp)
parser.add_argument('--normfile',
                    help="JSON file containing base normalizers",
                    default=None)

args = parser.parse_args()

header = " ".join([f"'{i}'" for i in sys.argv])

system = loos.createSystem(args.system_file)
all_trajs = []
for t in args.traj_files:
    traj = loos.pyloos.Trajectory(t, system,
                                  skip=args.skip,
                                  stride=args.stride)
    all_trajs.append(traj)

selection = loos.selectAtoms(system, args.selection_string)
if args.base_carbons:
    selection = loos.selectAtoms(selection, '!backbone && name =~"^C"')

residues = selection.splitByResidue()

scores = np.zeros((len(residues), len(residues)))

# default box in case the system isn't actually periodic
box = loos.GCoord(1000., 1000., 1000.)

total_frames = 0
for traj in all_trajs:
    for frame in traj:
        if frame.isPeriodic():
            box = frame.periodicBox()
        for i in range(len(residues)-1):
            for j in range(i+1, len(residues)):
                p = residues[i].stacking(residues[j], box, 5.0)
                scores[i, j] += p
                scores[j, i] += p
        total_frames += 1
scores /= total_frames

# Normalize, if they asked for it
if args.normfile:
    with open(args.normfile) as jsfile:
        stack_values = json.load(jsfile)

        resnames = []
        for r in residues:
            resnames.append(r[0].resname())

        for i in range(len(scores)-1):
            for j in range(i+1, len(scores)):
                key = resnames[i] + "-" + resnames[j]
                newval = scores[i, j] / stack_values[key]
                scores[i, j] = newval
                scores[j, i] = newval

np.savetxt(args.outfile_name, scores, header=header)
