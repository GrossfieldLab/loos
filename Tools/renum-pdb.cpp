/*
  renum-pdb

  Renumbers a PDB (though it could take any arbitrary model that LOOS
  understands...)
*/



/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo
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

typedef vector<AtomicGroup>  vGroup;


void show_help(void) {
  cout << "Usage- renum-pdb model selection resid-start atomid-start [selection resid-start atomid-start ...]\n";
  cout << "\n"
    "This tool renumbers both atomids and residue ids for a set of selections (independently).\n"
    "Note that atomids are important to LOOS when used with a trajectory.  If a model is\n"
    "renumbered, then it will no longer match the corresponding trajectory.\n";
  exit(0);
}



int main(int argc, char *argv[]) {

  // lightweight args processing...
  if (argc < 2 || (argc-2)%3)
    show_help();
  
  string name(argv[1]);
  if (name == "-h" || name == "--help")
    show_help();

  AtomicGroup model = createSystem(name);
  for (int i=2; i<argc; i += 3) {
    AtomicGroup subset = selectAtoms(model, argv[i]);
    vGroup residues = subset.splitByResidue();
    int resid = strtol(argv[i+1], (char **)NULL, 10);
    int atomid = strtol(argv[i+2], (char **)NULL, 10);

    vGroup::iterator vi;
    for (vi = residues.begin(); vi != residues.end(); ++vi, ++resid) {
      AtomicGroup::iterator gi;
      for (gi = (*vi).begin(); gi != (*vi).end(); ++gi, ++atomid) {
        (**gi).id(atomid);
        (**gi).resid(resid);
      }
    }
  }

  PDB pdb = PDB::fromAtomicGroup(model);
  cout << pdb;
}
