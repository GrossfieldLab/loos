#!/usr/bin/env python
"""
  axis_with_membrane computes the cosine between the first
  principal axis of a selection and the Z-axis (i.e. the putative membrane
  normal).  This is written out for a trajectory as a time series...
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
from math import *
import sys

## Command line arguments
model_name = sys.argv[1]
traj_name = sys.argv[2]
selection = sys.argv[3]

# Define the system first
mol = createSystem(model_name)

# Now define the trajectory
traj = createTrajectory(traj_name, mol)

# Select the atoms to be used for the principal axis calculation
backbone = selectAtoms(mol, selection)

# Iterate over all frames in the trajectory...
t = 0
while (traj.readFrame()):

    # Update the coordinates
    traj.updateGroupCoords(mol)
    
    # Compute the principal axes for the subset of atoms given above
    axes = backbone.principalAxes()

    # Print out time, and the dot-product between
    # The Z-axis (i.e. membrane normal) and the first principal axes
    # of the subset of atoms
    print t, "\t", axes[0].z()
    t += 1


