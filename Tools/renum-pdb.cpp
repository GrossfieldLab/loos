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




string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tRenumbers atoms and residues\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool renumbers sets of atoms and residues.  For each selection, the atomids and\n"
    "resids are incremented.  The rest of the model is left unchanged.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\trenum-pdb model.pdb 'all' 1 1 >renumbered.pdb\n"
    "This example renumbers everything, begin with resid 1 and atomid 1.\n"
    "\n"
    "\trenum-pdb model.pdb 'resid >= 100' 500 5000 >renumbered.pdb\n"
    "This example renumbers residues higher than 100, shifting them to begin with 500.\n"
    "The atomids for these residues are also renumbered, beginning with 5000.\n"
    "\n";

  return(msg);
}


int main(int argc, char *argv[]) {

  // lightweight args processing...
  if (argc < 2 || (argc-2)%3) {
    cerr << "Usage- renum-pdb model selection resid-start atomid-start [selection resid-start atomid-start ...]\n";
    cerr << fullHelpMessage();
    exit(-1);
  }
  
  string name(argv[1]);
  if (name == "-h" || name == "--help") {
    cerr << "Usage- renum-pdb model selection resid-start atomid-start [selection resid-start atomid-start ...]\n";
    cerr << fullHelpMessage();
    exit(-1);
  }
    

  AtomicGroup model = createSystem(name);
  for (int i=2; i<argc; i += 3) {
    AtomicGroup subset = selectAtoms(model, argv[i]);
    vGroup residues = subset.splitByResidue();
    int resid = strtol(argv[i+1], (char **)NULL, 10);
    int atomid = strtol(argv[i+2], (char **)NULL, 10);

    // Use AtomicGroup::renumber() to renumber the atomid's since this
    // will preserve the connectivity.
    subset.renumber(atomid);

    // Must manually update the resid's however...
    for (vGroup::iterator residue = residues.begin(); residue != residues.end(); ++residue) {
        for (AtomicGroup::iterator atom = residue->begin(); atom != residue->end(); ++atom)
            (*atom)->resid(resid);
        ++resid;
    }
    
  }

  PDB pdb = PDB::fromAtomicGroup(model);
  cout << pdb;
}
