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


pAtom findMatch(const pAtom& probe, const AtomicGroup& grp) {
  for (AtomicGroup::const_iterator i = grp.begin(); i != grp.end(); ++i)
    if ((*i)->name() == probe->name() && (*i)->id() == probe->id()
        && (*i)->resname() == probe->resname() && (*i)->resid() == probe->resid()
        && (*i)->segid() == probe->segid())
      return(*i);

  pAtom null;
  return(null);
}



int main(int argc, char *argv[]) {

  if (argc == 1) {
    cout << "Usage- side-nodes selection model [psf] >output.pdb\n";
    exit(0);
  }

  string hdr = invocationHeader(argc, argv);
  int k = 1;
  string selection(argv[k++]);
  AtomicGroup model = createSystem(argv[k++]);
  AtomicGroup subset = selectAtoms(model, selection);
  if (argc > k) {
    AtomicGroup structure = createSystem(argv[k++]);
    for (AtomicGroup::iterator i = subset.begin(); i != subset.end(); ++i) {
      pAtom match = findMatch(*i, structure);
      if (!match) {
        cerr << "ERROR- no match found for atom " << **i << endl;
        exit(-1);
      }
      
      if (!match->checkProperty(Atom::massbit)) {
        cerr << "ERROR- Atom has no mass: " << *match << endl;
        exit(-1);
      }
      
      (*i)->mass(match->mass());
    }
  }

  vector<AtomicGroup> residues = subset.splitByResidue();
  AtomicGroup cg_sites;
  int currid = model.maxId();

  for (vector<AtomicGroup>::iterator vi = residues.begin(); vi != residues.end(); ++vi) {
    // First, pick off the CA
    AtomicGroup CA = vi->select(AtomNameSelector("CA"));
    if (CA.empty()) {
      cerr << "Error- cannot find CA.\n" << *vi;
      exit(-10);
    }
    CA[0]->occupancy(CA[0]->mass());
    cg_sites += CA[0];

    AtomicGroup sidechain = vi->select(NotSelector(BackboneSelector()));
    if (sidechain.empty()) {
      cerr << "Warning- No sidechain atoms for:\n" << *vi;
      continue;
    }
    
    GCoord c = sidechain.centerOfMass();
    pAtom pa(new Atom(++currid, "CGS", c));
    pa->resid(CA[0]->resid());
    pa->resname(CA[0]->resname());
    pa->segid(CA[0]->segid());
    double m = sidechain.totalMass();
    pa->occupancy(m);

    cg_sites += pa;
  }

  PDB pdb = PDB::fromAtomicGroup(cg_sites);
  pdb.remarks().add(hdr);
  cout << pdb;
}
