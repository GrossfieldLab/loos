# !/usr/bin/env python3

"""
rmsds-align.py: Calcuate 2D RMSD using different alignment and rmsd selections.


The purpose of this tool is to calcuate pairwise RMSD for a selection after
aligning using a different selection. This script supports upto 2 trajectories.

Anees Mohammmed Keedakkatt Puthenpeedikakkal, 2022
"""

"""

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2013 Tod Romo, Grossfield Lab
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
import argparse
import loos
import loos.pyloos


def fullhelp():
    print("""
SYNOPSIS

Calculate a pair-wise RMSD for a trajectory (or two trajectories) with
different align and rmsd selections.


DESCRIPTION

This tool calculates the pairwise RMSD between each structures in a trajectory
after aligning using a selection which is different from the RMSD selection.
This tool can be used for two different trajectories also. This tool is similar
to rmsds, except rmsds uses same selection for alignment and pairwise
RMSD calculation.

The transformation matrix that aligns align_selection_2 is used to transform
align_selection_1, followed by RMSD calcuation using rmsd_selection_1 to
rmsd_selection_2. This is looped over the all the frames in the
trajectory (or for selected --range1 and --range2)

EXAMPLES

python3 rmsds-align.py model.pdb simulation.dcd > rmsd.asc

This example uses all alpha-carbons and every frame in the trajectory
(equivalent to rmsds).


python3 rmsds-align.py --model2 inactive.pdb --traj2 inactive.dcd active.pdb
active.dcd > rmsd.asc

This example uses all alpha-carbons and compares the "inactive" simulation
with the "active" one.


python3 rmsds-align.py  --align 'resid <= 20 && name== "CA"' --rmsd 'resid > 25
 && name =~ "^(C|CA|CB|O|N)$"' model.pdb simulation.dcd > rmsd.asc

This example uses only first 20 alpha carbons for alignment and pairwise RMSD
is calculated for all the backbone atoms with resid > 25.


python3 rmsds-align.py --align 'resid <= 20 && name== "CA"' --align2
'resid >=25 && resid <=44 && name == "CA"' --rmsd 'resid > 25 && resid < 50 &&
name =~ "^(C|CA|CB|O|N)$"' --rmsd2 'resid < 25 && name =~ "^(C|CA|CB|O|N)$"'
--model2 inactive.pdb --traj2 inactive.dcd active.pdb active.dcd > rmsd.asc

This complex example first aligns active’s first 20 alpha-Cs with inactive’s
alpha-Cs with resid >=25 and resid <=44. Pairwise RMSD is calculated between
active’s backbone (with resid’s in 25-49) and inactive’s backbone with
resid <25.


NOTES

When using two trajectories, the selections must match both in number of atoms
selected and in the sequence of atoms (i.e. the first atom in the align/rmsd
selection is matched with the first atom in the align2/rmsd2 selection.)

SEE ALSO
rmsd2ref, rmsds

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


header = " ".join(sys.argv)
parser = argparse.ArgumentParser(description="Calculate 2D RMSD using"
                " different alignment and rmsd selections.")

parser.add_argument('model', help="Model")

parser.add_argument('traj', help="Trajectory")

parser.add_argument('--model2',
                    help="Model for 2nd trajectory (default = model)",
                    default=None)

parser.add_argument('--traj2',
                    help="2nd trajectory (default = traj)",
                    default=None)

parser.add_argument('--align',
                    help="Align selection for 1st trajectory (default = CA)",
                    default='name == "CA"', type=str)

parser.add_argument('--align2',
                    help="Align selection 2nd trajectory (default = align)",
                    default=None, type=str)

parser.add_argument('--rmsd',
                    help="RMSD selection for 1st trajectory"
                    "(default = (C|O|N|CA))",
                    default='name =~ "^(C|O|N|CA)$"', type=str)

parser.add_argument('--rmsd2',
                    help="RMSD selection for 2nd trajectory (default = rmsd)",
                    default=None, type=str)

parser.add_argument('--precision',
                    help="Write out matrix coefficients with this many digits "
                    "(default = 2)", default=2, type=int)

parser.add_argument('--range1',
                    help="Range of frames to use from the first trajectory "
                    "[FORMAT - skip:stride:stop, stop is optional] "
                    "(default = 0:1:last)", default='0:1', type=str)

parser.add_argument('--range2', help="Range of frames to use from the second "
                    "trajectory [FORMAT - skip:stride:stop, stop is optional] "
                    "(default = range1)",
                    default=None, type=str)

parser.add_argument('--fullhelp',
                    help="Print detailed description of all options",
                    action=FullHelp)

args = parser.parse_args()

# In case of only one trajectory is given as input
if args.rmsd2 is None:
    args.rmsd2 = args.rmsd

if args.align2 is None:
    args.align2 = args.align

if args.model2 is None:
    args.model2 = args.model

if args.traj2 is None:
    args.traj2 = args.traj

if args.range2 is None:
    args.range2 = args.range1


# Organizing range of frames to consider in trajectories
range_1 = args.range1.split(":")
range_2 = args.range2.split(":")

# In case only skip and stride is given as input for traj
if len(range_1) == 2:
    select_1 = range(int(range_1[0]), len(args.traj), int(range_1[1]))
# In case only skip, stride and stop are given
elif len(range_1) == 3:
    select_1 = range(int(range_1[0]), int(range_1[2]), int(range_1[1]))
else:
    print("Error in range")
    sys.exit(0)


# In case only skip and stride is given as input for traj2
if len(range_1) == 2:
    select_2 = range(int(range_2[0]), len(args.traj2), int(range_2[1]))
# In case only skip, stride and stop are given
elif len(range_1) == 3:
    select_2 = range(int(range_2[0]), int(range_2[2]), int(range_2[1]))
else:
    print("Error in range")
    sys.exit(0)

# Model and Trajectory
model = loos.createSystem(args.model)
traj = loos.pyloos.Trajectory(args.traj, model, iterator=select_1)
model2 = loos.createSystem(args.model2)
traj2 = loos.pyloos.Trajectory(args.traj2, model2, iterator=select_2)

# Align and RMSD selection definitions
align_selection_1 = loos.selectAtoms(model, args.align)
align_selection_2 = loos.selectAtoms(model2, args.align2)

rmsd_selection_1 = loos.selectAtoms(model, args.rmsd)
rmsd_selection_2 = loos.selectAtoms(model2, args.rmsd2)

# Print input
print("# Align1 - ", args.align)
print("# Align2 - ", args.align2)
print("# RMSD1 - ", args.rmsd)
print("# RMSD2 - ", args.rmsd2)
print("# Traj-range1 - ", args.range1)
print("# Traj-range2 - ", args.range2)


for frame in traj:
    ref_align = align_selection_1.copy()
    ref_target = rmsd_selection_1.copy()

    for frame in traj2:
        # Find transformation matrix that aligns align_selection_2
        # onto ref_align(align_selection_1)
        trans_matrix = ref_align.alignOnto(align_selection_2)
        # Tranform matrix
        transform = loos.XForm(trans_matrix)
        # Apply the tranform to ref_target(rmsd_selection_1)
        ref_target.applyTransform(transform)
        # Calculate the RMSD between rmsd_selection_2 and
        # ref_target(rmsd_selection_1)
        rmsd_value = rmsd_selection_2.rmsd(ref_target)
        # Print RMSD
        print(round(rmsd_value, args.precision), "", end='')

    print("")
