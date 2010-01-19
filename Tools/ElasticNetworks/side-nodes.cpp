/*
  side-nodes

  (c) 2010 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry


  Usage:
    side-nodes selection model >output.pdb

*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009-2010 Tod D. Romo
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


using namespace std;
using namespace loos;

int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);
  int k = 1;
  string selection(argv[k++]);
  AtomicGroup model = createSystem(argv[k++]);
  AtomicGroup subset = selectAtoms(model, selection);

  vector<AtomicGroup> residues = subset.splitByResidue();
  AtomicGroup cg_sites;
  int currid = model.maxId();

  for (vector<AtomicGroup>::iterator vi = residues.begin(); vi != residues.end(); ++vi) {
    // First, pick off the CA
    AtomicGroup CA = (*vi).select(AtomNameSelector("CA"));
    if (CA.empty()) {
      cerr << "Error- cannot find CA.\n" << *vi;
      exit(-10);
    }
    cg_sites += CA[0];

    AtomicGroup sidechain = (*vi).select(NotSelector(BackboneSelector()));
    if (sidechain.empty()) {
      cerr << "Error- No sidechain atoms for:\n" << *vi;
      exit(-10);
    }
    
    GCoord c = sidechain.centerOfMass();
    pAtom pa(new Atom(++currid, "CGS", c));
    pa->resid(CA[0]->resid());
    pa->resname(CA[0]->resname());
    pa->segid(CA[0]->segid());

    cg_sites += pa;
  }

  PDB pdb = PDB::fromAtomicGroup(cg_sites);
  pdb.remarks().add(hdr);
  cout << pdb;
}
