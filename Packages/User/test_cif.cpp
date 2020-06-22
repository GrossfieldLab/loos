/*
  simple_model_calc.cpp

  (c) 2011 Tod D. Romo, Grossfield Lab
           Department of Biochemistry
           University of Rochster School of Medicine and Dentistry


  C++ template for writing a tool that performs a calculation on a model with minimal
  command-line options
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011 Tod D. Romo
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

  // How the tool was invoked, for logging purposes...
  string header = invocationHeader(argc, argv);

  // ***EDIT***
  // Verify correct number of command-line arguments
  if (argc != 3) {
    cerr << "Usage- simple_model_calc model-name selection\n";
  }

  // Handle command-line arguments
  int arg_index = 1;

  // Read in the model
  AtomicGroup model = createSystem(argv[arg_index++]);

  PDB pdb = PDB::fromAtomicGroup(model);
  MMCIF cif = MMCIF::fromAtomicGroup(model);

  cout << cif << endl;

}
