#!/usr/bin/env python3

import sys
import loos
import LipidLibrary
import random


protein_pdb = sys.argv[1]
num_ligands = int(sys.argv[2])
box_size = float(sys.argv[3])
library_location = sys.argv[4]
output_pdb = sys.argv[5]

box = loos.GCoord(box_size, box_size, box_size)
half_box = 0.5 * box_size

protein = loos.createSystem(protein_pdb)
protein.periodicBox(box)

library = LipidLibrary.LipidLibrary(library_location)

# Might want to make this optional
protein.centerAtOrigin()

accepts = 0
trials = 0
while accepts < num_ligands:
    trials += 1

    new_molecule = library.pick_structure()
    new_molecule.centerAtOrigin()

    # translate by random within box_size
    x_trans = random.uniform(-half_box, half_box)
    y_trans = random.uniform(-half_box, half_box)
    z_trans = random.uniform(-half_box, half_box)
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
      accepts/trials, " %")

# protein.renumber()
pdb = loos.PDB.fromAtomicGroup(protein)

with open(output_pdb, "w") as outfile:
    outfile.write(str(pdb))
