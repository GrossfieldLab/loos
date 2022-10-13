#!/usr/bin/env python3
"""
Insert a number of small molecules into a simulation of a protein or RNA.
This can be viewed as preparatory to running solvate.py or omg.py

Alan Grossfield,  University of Rochester Medical Center, 2019
"""

"""

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2019 Tod Romo, Grossfield Lab
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
import loos.OptimalMembraneGenerator
from loos.OptimalMembraneGenerator import LipidLibrary
import random
import argparse

def main():
    parser = argparse.ArgumentParser(description="Add small molecules to a system")

    parser.add_argument('num_ligands',
                        help="Number of small molecules",
                        type=int)
    parser.add_argument('box_size',
                        help="Dimension of cubic box, in Angstroms",
                        type=float)
    parser.add_argument('library_location',
                        help="Directory with library of small molecules")
    parser.add_argument('output_pdb',
                        help="Name of output PDB file")
    parser.add_argument('--protein',
                        help="File containing coordinates of the protein to be surrounded by small molecules")
    parser.add_argument('--no_center',
                        help="If set, do not center the protein at the origin",
                        action="store_true")
    parser.add_argument('--zbox',
                        help="Z dimension of system, if not cubic",
                        type=float)
    parser.add_argument('--z_exclude',
                        help="Don't let molecules be placed in +/- z_exclude",
                        type=float)
    args = parser.parse_args()


    box = loos.GCoord(args.box_size, args.box_size, args.box_size)
    half_box = 0.5 * args.box_size
    half_z = half_box

    if args.zbox:
        box.z(args.zbox)
        half_z = 0.5 * args.zbox

    if args.z_exclude:
        args.z_exclude = abs(args.z_exclude)
    else:
        args.z_exclude = 0.0


    if args.protein:
        protein = loos.createSystem(args.protein)
    else:
        protein = loos.AtomicGroup()

    protein.periodicBox(box)

    library = LipidLibrary.LipidLibrary(args.library_location)

    if not args.no_center:
        protein.centerAtOrigin()

    accepts = 0
    trials = 0
    while accepts < args.num_ligands:
        trials += 1

        new_molecule = library.pick_structure()
        new_molecule.centerAtOrigin()

        # translate by random within box_size
        x_trans = random.uniform(-half_box, half_box)
        y_trans = random.uniform(-half_box, half_box)
        z_trans = random.uniform(-half_z, half_z)

        if abs(z_trans) < args.z_exclude:
            continue

        trans = loos.GCoord(x_trans, y_trans, z_trans)

        new_molecule.translate(trans)

        # bump check against existing system
        if not protein.contactWith(6., new_molecule, box):
            # fix the residue number
            accepts += 1
            for i in range(len(new_molecule)):
                new_molecule[i].resid(accepts)
            protein.append(new_molecule)

    print("Placed ", accepts, " molecules in ", trials, " trials: ",
        accepts/trials * 100, " %")

    protein.renumber()
    protein.clearBonds()
    pdb = loos.PDB.fromAtomicGroup(protein)

    with open(args.output_pdb, "w") as outfile:
        outfile.write(str(pdb))

if __name__ == "__main__":
    main()
    