/*
  model-select.cpp

  Takes a model (PDB, PSF, etc) and a selection string.  Parses the selection, then
  applies it to the PDB and writes the output to stdout.  This tool is
  used maily for checking your selection strings to make sure you're
  actually selecting what you intend to select...
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

using namespace std;
using namespace loos;


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tRaw dump of a model subset in LOOS\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool is useful for diagnosing problems with selections and how\n"
    "LOOS reads model files.  It will write out a pseudo-XML representation\n"
    "of the information it has stored about the selected subset.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tmodel-select all model.pdb >model.xml\n"
    "This example writes out ALL atoms\n"
    "\n"
    "\tmodel-select 'name == \"CA\"' model.pdb >model-ca.xml\n"
    "This example only writes out alpha-carbons.\n"
    "\n";

  return(msg);
}



int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);

  if (argc != 3) {
    cerr << "Usage- " << argv[0] << " selection-string model\n";
    cerr << fullHelpMessage();
    exit(-1);
  }

  AtomicGroup model = createSystem(argv[2]);
  AtomicGroup subset = selectAtoms(model, argv[1]);

  cerr << "You selected " << subset.size() << " atoms out of " << model.size() << endl;
  
  cout << "<!-- " << header << " -->\n";
  cout << subset << endl;
}
