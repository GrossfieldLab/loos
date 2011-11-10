/*
  water-count.cpp

   usage:
     water-count prefix >output.asc

   Description:
     Counts the number of waters that are inside the protein at each
     timepoint and writes this out as a vector.  In other words, it
     just sums the rows for each column of the water matrix.
*/



/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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


using namespace std;
using namespace loos;


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  if (argc != 2) {
    cerr << "Usage - water-count prefix >output.asc\n";
    exit(-1);
  }

  string prefix(argv[1]);

  Math::Matrix<double> V;
  readAsciiMatrix(prefix + ".vol", V);

  Math::Matrix<int> M;
  readAsciiMatrix(prefix + ".asc", M);
  uint m = M.rows();
  uint n = M.cols();

  if (V.rows() != M.cols()) {
    cerr << "ERROR - mismatch in volume and water matrix\n";
    exit(-10);
  }

  cout << "# " << hdr << endl;
  cout << "# frame\tcount\tvolume\n";
  for (uint i=0; i<n; ++i) {
    long cnt = 0;
    for (uint j=0; j<m; ++j)
      cnt += M(j,i);
    cout << i << "\t" << cnt << "\t" << V(i,0) << endl;
  }
}
