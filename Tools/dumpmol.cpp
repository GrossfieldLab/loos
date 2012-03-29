/*
  dumpmol.cpp

  Reads in a system model (PDB, PSF, etc) and writes out the AtomicGroup
  base representation.  This tool is useful for checking to make sure that
  the input model is being correctly parsed by LOOS.
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
  std::cerr << "WARNING- this tool is deprecated and will be removed in a future release\n";
  std::cerr << "         Use the model-select tool instead.\n";

  if (argc != 2) {
    std::cerr << "Usage - dumpmol model\n";
    exit(-1);
  }

  AtomicGroup model = createSystem(argv[1]);
  std::cout << model << std::endl;
}
