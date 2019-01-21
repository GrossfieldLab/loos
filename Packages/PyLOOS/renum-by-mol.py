#!/bin/env python
"""
  Renumber a model so that the resids for each molecule start at 1
"""
"""

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2016 Tod Romo, Grossfield Lab
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


if len(sys.argv) == 1 or sys.argv[1] == '-h':
    print('Usage- renum-by-mol.py model [selection]')
    sys.exit(0)
    

model = loos.createSystem(sys.argv[1])

if len(sys.argv) > 2:
    model = loos.selectAtoms(model, sys.argv[2])
    model.pruneBonds()

if model.hasBonds():
    molecules = model.splitByMolecule()
else:
    molecules = model.splitByUniqSegid()

for molecule in molecules:
    residues = molecule.splitByResidue()
    resid = 1
    for residue in residues:
        for atom in residue:
            atom.resid(resid)
        resid += 1

pdb = loos.PDB.fromAtomicGroup(model)
print(str(pdb), end=' ')
