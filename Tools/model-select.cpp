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

using namespace loos;

int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);

  if (argc != 3) {
    cerr << "Usage- " << argv[0] << " <selection string> <pdb file>\n";
    exit(-1);
  }

  AtomicGroup model = createSystem(argv[2]);
  AtomicGroup subset = selectAtoms(model, argv[1]);

  cerr << "You selected " << subset.size() << " atoms out of " << model.size() << endl;
  
  cout << "<!-- " << header << " -->\n";
  cout << subset << endl;
}
