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




#include <pdb.hpp>
#include <Parser.hpp>
#include <Selectors.hpp>


int main(int argc, char *argv[]) {
  
  PDB p(argv[1]);

  cout << "Read in " << p.size() << " atoms from " << argv[1] << endl;

  Parser parsed(argv[2]);
  KernelSelector parsed_sel(parsed.kernel());
  AtomicGroup  usel = p.select(parsed_sel);
  
  cout << "There are " << usel.size() << " atoms in the selection.\n";
  vector<AtomicGroup> grps = usel.splitByUniqueSegid();

  cout << "There are " << grps.size() << " groups selected.\n";

  cout << grps[0] << endl;

}

