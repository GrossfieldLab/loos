#!/usr/bin/env python3

"""

Given a model file that has connectivity and a file that has good metadata
write a new pdb file that has the metadata and connectivity information.

Example usage case would be combining a Tinker XYZ file and equivalent PDB file
to produce a PDB file with bond information.

THIS PROGRAM ASSUMES THE CONTENTS OF THE TWO MODELS ARE IDENTICAL: the same
atoms in the same order.

"""

"""

    This file is part of LOOS.

    LOOS (Lightweight Object-Oriented Structure library)
    Copyright (c) 2022 Alan Grossfield
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
import loos.pyloos
import argparse
import sys

parser = argparse.ArgumentParser(description="Add bond info to a PDB")
parser.add_argument('main_model_file',
                    help="Model file without bond info")
parser.add_argument('model_file_with_bonds',
                    help="Model file with bond info")
parser.add_argument('new_pdb',
                    help="Output PDB file, with bonds")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)
args = parser.parse_args()

main = loos.createSystem(args.main_model_file)
with_bonds = loos.createSystem(args.model_file_with_bonds)

# do some sanity checking
if len(main) != len(with_bonds):
    print("The two models must have the same number of atoms, in the same order")
    print(f"Main file: {len(main)} atoms")
    print(f"Bonds file: {len(with_bonds)} atoms")
    sys.exit(-1)

if not with_bonds.hasBonds():
    print(f"Second file ({args.model_file_with_bonds}) must have bond information")
    sys.exit(-1)

for i in range(len(main)):
    m = main[i]
    w = with_bonds[i]

    #if m.name() != w.name():
    #    print(f"Names don't match for atom {i}: {m.name()} {w.name()}")
    #    sys.exit(-1)

    for bonded in w.getBonds():
        m.addBond(bonded)

new = loos.PDB.fromAtomicGroup(main)
with open(args.new_pdb, "w") as out_pdb:
    out_pdb.write(str(new))
