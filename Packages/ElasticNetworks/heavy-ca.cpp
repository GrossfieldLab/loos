/*
  heavy-ca

  Given a PDB where masses are stored in the occupancy field, reduce
  the structure to CA's only where the mass of the CA is the sum of
  the mass of all atoms in the corresponding residue...

  Usage- heavy-ca selection model >output
      
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009 Tod D. Romo
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



string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "Store whole residue mass in CA occupancy field\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "Given a PDB where masses are stored in the occupancy field, reduce\n"
    "the structure to CA's only where the mass of the CA is the sum of\n"
    "the mass of all atoms in the corresponding residue.\n"
    "\n"
    "Note: The selection string in this tool is used to decide which\n"
    "      residues to sum the mass of.  So 'name==\"CA\"' will return\n"
    "      a mass of 12 in the occupancy field.  \n"
    "\n"
    //
    "EXAMPLES\n"
    "\n"
    "heavy-ca 'segid==\"PROT\"' model.pdb > newmodel.pdb'\n"
    "\tMasses in the occupancy field of model.pdb are \n"
    "\tsummed over each residue in segid PROT and placed\n"
    "\ton the CA in newmodel.pdb.\n"
    "\n"
    "heavy-ca 'segid==\"PROT\" && !(hydrogen)' model.pdb > newmodel.pdb'\n"
    "\tSame as above, but hydrogen atoms are excluded \n"
    "\tfrom the summation.\n"
    "\n"
    "SEE ALSO\n"
    "\n"
    "Packages/ElasticNetworks/psf-masses - \n"
    "\tThis tool will take the masses from a PSF file\n"
    "\tand place them in the occupancy field of a PDB\n"
    "\n"
    "\n";

  return(msg);
}



int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage- heavy-ca selection pdb >output\n";
    cerr << fullHelpMessage();
    exit(0);
  }


  string hdr = invocationHeader(argc, argv);

  string selection(argv[1]);
  AtomicGroup structure = createSystem(argv[2]);

  AtomicGroup model = selectAtoms(structure, selection);

  vector<AtomicGroup> residues = model.splitByResidue();
  AtomicGroup heavy;

  for (vector<AtomicGroup>::iterator j = residues.begin(); j != residues.end(); ++j) {
    AtomicGroup grp = (*j).select(AtomNameSelector("CA"));
    if (grp.empty()) {
      cerr << "ERROR- could not find a CA in the following residue:\n";
      cerr << *j;
      exit(-1);
    }

    pAtom ca = grp[0];
    double mass = 0.0;
    for (AtomicGroup::iterator i = (*j).begin(); i != (*j).end(); ++i)
      mass += (*i)->occupancy();

    ca->occupancy(mass);
    heavy.append(ca);
  }


  PDB outpdb = PDB::fromAtomicGroup(heavy);
  outpdb.remarks().add(hdr);
  cout << outpdb;
}
