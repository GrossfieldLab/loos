#!/usr/bin/env python3
"""
Print a bunch of useful information about the system to stdout

Alan Grossfield, 2021
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
import loos

filename = sys.argv[1]

system = loos.createSystem(filename)

print(f"Num atoms = {len(system)}")
print(f"Has Coords = {system.hasCoords()}")

hasCharge = system[0].checkProperty(loos.Atom.chargebit)
print(f"Has Charges = {hasCharge}")
if hasCharge:
    print(f"Total charge = {system.totalCharge():.4f}")

print(f"Has Masses = {system[0].checkProperty(loos.Atom.massbit)}")

print(f"Has Bonds = {system.hasBonds()}")
print(f"Is Periodic = {system.isPeriodic()}")
if system.isPeriodic():
    print(f"Box = {system.periodicBox()}")

if system.hasCoords():
    print(f"Centroid = {system.centroid()}")
