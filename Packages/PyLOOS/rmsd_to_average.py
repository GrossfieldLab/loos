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
from loos import *
import sys

## Command line arguments
model_name = sys.argv[1]
traj_name = sys.argv[2]
align_with = sys.argv[3]
rmsd_with = sys.argv[4]
skip = int(sys.argv[5])

# Create the model & read in the trajectory
model = createSystem(model_name)
traj = createTrajectory(traj_name, model)

align_subset = selectAtoms(model, align_with)
rmsd_subset = selectAtoms(model, rmsd_with)

print "# Alignment ubset has %d atoms." % (len(align_subset))

patraj = PyAlignedTraj(traj, rmsd_subset, skip = skip, alignwith = align_subset)
average = patraj.averageStructure()

t = 0
avg_rmsd = 0
for structure in patraj:
   rmsd = average.rmsd(structure)
   avg_rmsd = avg_rmsd + rmsd
   print "%d\t%f" % (t, rmsd)
   t = t + 1

print "# Average rmsd = %f" % (avg_rmsd/t)



