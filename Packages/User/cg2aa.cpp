/*
  cg2aa.cpp


  Converts a LOOS-supported format to a PDB (so long as coordinates
  are present)

  Usage:

    cg2aa structure-file mapping-database >output.pdb

*/




/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
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
*/

#include <loos.hpp>
#include <XForm.hpp>

using namespace std;
using namespace loos;


int main(int argc, char *argv[]) {
  
  // Read in structure
  AtomicGroup model = createSystem(argv[1]);

  // Read in mapping scheme
  AtomicGroup cgmap = createSystem(argv[2]);
  AtomicGroup aamap = createSystem(argv[3]);

  // Iterate over residues of given input...
  vector<AtomicGroup> inputModel;
  inputModel = model.splitByResidue();

  AtomicGroup output;

  for (vector<AtomicGroup>::iterator m =inputModel.begin(); m!=inputModel.end(); ++m)
      {

      string resn = (m->getAtom(0))->resname();
      string match = string("resname == \"" + resn + "\"");
      
      AtomicGroup cgMatch = selectAtoms(cgmap, match);
      GMatrix X = cgMatch.superposition(*m);
      XForm Xa;
      Xa.load(X);
      AtomicGroup aaMatch = selectAtoms(aamap, match);
      AtomicGroup temp = aaMatch.copy();
      temp.applyTransform(Xa);
      output.append(temp);

      }

// superposition and transformation
// GMatrix X = cgGroup.superposition(cginputGroup);
// XForm Xa = XForm(X);
// copy the aa version of the group (deep copy), call it temp
// AtomicGroup temp = aaGroup.applyTransform(Xa);
// append temp to running AtomicGroup

  PDB pdb = PDB::fromAtomicGroup(output);
  pdb.remarks().add(invocationHeader(argc, argv));

  cout << pdb;
}
