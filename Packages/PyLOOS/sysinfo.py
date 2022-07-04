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
import argparse

parser = argparse.ArgumentParser(description="Print system information")
parser.add_argument('system_file', help="File describing the system")

# option to output json
parser.add_argument('--json', '-j',
                    default=None,
                    help="Specify a file for json output")
args = parser.parse_args()

system = loos.createSystem(args.system_file)

num_atoms = len(system)
has_coords = system.hasCoords()
has_charge = system[0].checkProperty(loos.Atom.chargebit)
if has_charge:
    total_charge = system.totalCharge()
else:
    total_charge = 0.0
has_mass = system[0].checkProperty(loos.Atom.massbit)
has_bonds = system.hasBonds()
is_periodic = system.isPeriodic()
box = system.periodicBox()
centroid = system.centroid()

if args.json:
    import json
    d = {}
    d["system_file"] = args.system_file
    d["num_atoms"] = num_atoms
    d["has_coords"] = has_coords
    d["has_charge"] = has_charge
    d["total_charge"] = total_charge
    d["has_mass"] = has_mass
    d["has_bonds"] = has_bonds
    d["is_periodic"] = is_periodic
    d["box"] = str(box)
    d["centroid"] = str(centroid)

    with open(args.json, "w") as json_file:
        json.dump(d, json_file, sort_keys=True)
else:
    print(f"System file = {args.system_file}")
    print(f"Num atoms = {num_atoms}")
    print(f"Has Coords = {has_coords}")
    print(f"Has Charges = {has_charge}")
    if has_charge:
        print(f"Total charge = {total_charge:.4f}")
    print(f"Has Masses = {has_mass}")
    print(f"Has Bonds = {has_bonds}")
    print(f"Is Periodic = {is_periodic}")
    if system.isPeriodic():
        print(f"Box = {box}")

    if system.hasCoords():
        print(f"Centroid = {centroid}")
