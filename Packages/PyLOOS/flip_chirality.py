#!/usr/bin/env python
"""
  flip_chirality is a largely pointless tool -- it reads in a structure,
  reverses the sign of the x coordinate, and prints it as a PDB to stdout.
"""
"""

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2012 Alan Grossfield, Grossfield Lab
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

import loos
import sys

filename = sys.argv[1]

system = loos.createSystem(filename)

system.centerAtOrigin()
# loop over atoms and reverse the sign of the x coordinate
for i in range(system.size()):
    # 2 different ways to assign to coordinates
    # SWIG won't let us do: coords.x() = new_x
    # so we can do this:
    #system[i].coords().x(-system[i].coords().x())
    # or this:
    system[i].coords()[0] *= -1.0

system.centerAtOrigin()

# convert to PDB and print out
pdb = loos.PDB_fromAtomicGroup(system)

# now rotate 180 degrees around the z axis
z_axis = loos.GCoord(0,0,1)
pdb.rotate(z_axis, 180)

# print the pdb out
print pdb



