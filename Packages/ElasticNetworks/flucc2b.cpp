/*
  flucc2b


   Assign fluctuations to a PDB...
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008 Tod D. Romo
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
#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>

#include <cassert>

using namespace std;
using namespace loos;

const double kB = 1.3606504e-23;   // Boltzmann

typedef Math::Matrix<double, Math::ColMajor> Matrix;

void show_help(void) {
  cerr << "Usage- flucc2b [selection] model-name pseudo-inverse scaling >output.pdb\n";
  exit(0);
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  cerr << "Warning- this tool is now deprecated.  Please use eigenflucc instead.\n";

  if (argc < 5 || argc > 6)
    show_help();

  int idx = 1;

  string selection("name == 'CA'");
  if (argc == 6)
    selection = string(argv[idx++]);

  string model_name(argv[idx++]);
  string pseudo_name(argv[idx++]);
  double scale = strtod(argv[idx++], 0);

  AtomicGroup model = createSystem(model_name);
  AtomicGroup subset = selectAtoms(model, selection);

  Matrix G;
  readAsciiMatrix(pseudo_name, G);
  int m = G.rows();
  int n = G.cols();

  if (m != n) {
      cerr << "ERROR - pseudoinverse matrix is non-square" << endl;
      exit(-1);
  }

  if (m == subset.size()) {

    // Assign B-factor (see Atilgan et al, Biophysical J. 2001 80:505-515, eq 8
    for (int i=0; i<subset.size(); ++i) {
      double d = scale * G(i,i);
      subset[i]->bfactor(d);
    }

  } else if (m == 3 * subset.size()) {  // This came from ANM...

    for (int i=0; i<subset.size(); ++i) {
      int j = i * 3;
      double d = scale * G(j,j) * G(j+1, j+1) * G(j+2, j+2);
      subset[i]->bfactor(d);
    }

  }
  
  // Up-cast AtomicGroup to a PDB
  PDB output = PDB::fromAtomicGroup(model);
  output.remarks().add(hdr);
  cout << output;

}
