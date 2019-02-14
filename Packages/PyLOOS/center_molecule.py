#!/usr/bin/env python3
"""
  Center an entire model based on a selection.
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

# Command line arguments
model_name = sys.argv[1]
selection = sys.argv[2]

# Create the model
model = createSystem(model_name)

# Select the requested atoms
subset = selectAtoms(model, selection)

# Compute centroid the old-fashioned way...
# Obviously, it is faster to call subset.centroid(),
# but we're demonstrating how else one could do it.
center = GCoord(0,0,0)
for atom in subset:
    center = center + atom.coords()
center = center / len(subset)

for atom in subset:
    atom.coords(atom.coords() - center)

# Convert to a PDB
pdb = PDB.fromAtomicGroup(model)

# Add a REMARK to the PDB
pdb.remarks().add("Structure centered using '" + selection + "'")

# Print it to stdout...
print(pdb)
