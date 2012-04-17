#  This file is part of LOOS.
#
#  LOOS (Lightweight Object-Oriented Structure library)
#  Copyright (c) 2012 Tod D. Romo, Grossfield Lab
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



# Import LOOS 
from loos import *
from math import *

# Define the system first
mol = createSystem("lfb.psf")

# Now define the trajectory
traj = createTrajectory("lfb.dcd", mol)

# Select the backbone atoms coming from segment PE1
backbone = selectAtoms(mol, "name =~ '^(C|O|N|CA)$' && segid == 'PE1'")

# Iterate over all frames in the trajectory...
t = 0
while (traj.readFrame()):

    # Update the coordinates
    traj.updateGroupCoords(mol)
    
    # Compute the principal axes for the subset of atoms given above
    axes = backbone.principalAxes()

    # Print out time, and the cosine of the dot-product between
    # The Z-axes (i.e. membrane normal) and the first principal axes
    # of the subset of atoms
    print t, "\t", cos(axes[0].z())
    t += 1


