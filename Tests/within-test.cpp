/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2012, Tod D. Romo
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

  AtomicGroup model = createSystem("./Tests/test.pdb");
  
  AtomicGroup lipid1 = selectAtoms(model, "segid == 'L1'");
  AtomicGroup palm = selectAtoms(lipid1, "resname == 'PALM'");
  AtomicGroup oleo = selectAtoms(lipid1, "resname == 'OLEO'");
  AtomicGroup lipid2 = selectAtoms(model, "segid == 'L2'");
  AtomicGroup protein = selectAtoms(model, "segid == 'A'");

  AtomicGroup result = palm.within(2.5, oleo);
  cout << "===> PALM -> OLEO\n";
  cout << "*** 2.5 angstroms\n";
  cout << result << endl;

  result = palm.within(10, oleo);
  cout << "*** 10 angstroms\n";
  cout << result << endl;
  
  result = lipid1.within(10, lipid2);
  cout << "==> lipid1 -> lipid2\n";
  cout << "*** 10 angstroms\n";
  cout << result << endl;

  result = lipid1.within(50, lipid2);
  cout << "*** 50 angstroms\n";
  cout << result << endl;
  
  result = lipid1.within(20, protein);
  cout << "==> lipid1 -> protein without box\n";
  cout << "*** 20 angstroms\n";
  cout << result << endl;


  result = lipid1.within(20, protein, GCoord(20,30,40));
  cout << "==> lipid1 -> protein with box\n";
  cout << "*** 20 angstroms\n";
  cout << result << endl;
    
}
