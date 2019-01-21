#!/usr/bin/env python
"""
  Compute an optimal average structure using an iterative alignment method, then
  output the RMSD from each frame and the average structure
"""
"""

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2012 Tod Romo, Grossfield Lab
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



# Import LOOS 
import loos
import loos.pyloos
import argparse
import sys


parser = argparse.ArgumentParser()
parser.add_argument('model', help='Structure to use')
parser.add_argument('traj', help='Trajectory')
parser.add_argument('align', help='Selection to use for aligning')
parser.add_argument('--skip', help='Skip this amount from the start of the trajectory', type=int, default=0)
parser.add_argument('--stride', help='Step through the trajectory by this many frames', type=int, default=1)
parser.add_argument('--rmsd', help='Use this selection for computing the RMSD', default='all')
args = parser.parse_args()

# Create the model & read in the trajectory
model = loos.createSystem(args.model)
traj = loos.pyloos.Trajectory(args.traj, model, skip = args.skip, stride = args.stride)

patraj = loos.pyloos.AlignedVirtualTrajectory(traj, alignwith = args.align)
patraj.setSubset(args.rmsd)
average = loos.pyloos.averageStructure(patraj)


avg_rmsd = 0
for structure in patraj:
   rmsd = average.rmsd(structure)
   avg_rmsd = avg_rmsd + rmsd
   print("%d\t%f" % (patraj.index(), rmsd))

print("# Average rmsd = %f" % (avg_rmsd/len(patraj)))



