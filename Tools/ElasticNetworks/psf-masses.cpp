/*
  psf-masses

  (c) 2009 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

  Takes masses from a PSF file and places them into the occupancy field of a PDB.
      
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008,2009 Tod D. Romo
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


using namespace loos;
using namespace std;


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
  if (argc != 3) {
    cerr << "Usage- psf-masses model.psf model.pdb >newmodel.pdb\n";
    exit(0);
  }

  string hdr = invocationHeader(argc, argv);

  AtomicGroup source = createSystem(argv[1]);
  AtomicGroup target = createSystem(argv[2]);

  for (AtomicGroup::iterator i = target.begin(); i != target.end(); ++i) {
    pAtom match = findMatch(*i, source);
    if (!match) {
      cerr << "ERROR- no match found for atom " << **i << endl;
      exit(-1);
    }

    if (!match->checkProperty(Atom::massbit)) {
      cerr << "ERROR- Atom has no mass: " << *match << endl;
      exit(-1);
    }

    (*i)->occupancy(match->mass());
  }


  PDB pdb = PDB::fromAtomicGroup(target);
  pdb.remarks().add(hdr);
  cout << pdb;
}
