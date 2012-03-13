/*
  gridslice.cpp

  Takes a double grid and extracts a plane from it as a matrix...
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
#include <boost/tuple/tuple.hpp>
#include <DensityGrid.hpp>

using namespace std;
using namespace loos;
using namespace loos::DensityTools;

typedef Math::Matrix<double, Math::RowMajor> Matrix;

void invalidIndex(int i) {
  cerr << "ERROR - invalid plane index " << i << endl;
  exit(-1);
}



int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage- gridslice [i|j|k] index <grid >matrix\n";
    cerr <<
      "\n"
      "Gridslice extracts a slice of the grid and writes it out\n"
      "as a Matlab/Octave/Gnuplot compatible ASCII matrix.\n"
      "The first option (i, j, or k) determines the orientation\n"
      "of the slice.  The index represents the coordinate in the\n"
      "direction.  For example, \"k 20\" means extract the plane\n"
      "when k=20 (an i,j-plane).  Using "i 13" means extract the\n"
      "plane when i=13 (a j,k-plane).\n";
    exit(-1);
  }

  string hdr = invocationHeader(argc, argv);

  string plane(argv[1]);
  int idx = atoi(argv[2]);

  DensityGrid<double> grid;
  cin >> grid;
  DensityGridpoint dims = grid.gridDims();
  cerr << boost::format("Grid dimensions are %d x %d x %d (i x j x k)\n") % dims[0] % dims[1] % dims[2];
  if (plane == "k") {

    if (idx > dims[2])
      invalidIndex(idx);
    Matrix M(dims[1]+1, dims[0]+1);
    for (int j=0; j<dims[1]; ++j)
      for (int i=0; i<dims[0]; ++i)
        M(j,i) = grid(idx,j,i);

    writeAsciiMatrix(cout, M, hdr);

  } else if (plane == "j") {

    if (idx > dims[1])
      invalidIndex(idx);
    Matrix M(dims[2]+1, dims[0]+1);
    for (int k=0; k<dims[2]; ++k)
      for (int i=0; i<dims[0]; ++i)
        M(k,i) = grid(k,idx,i);
    
    writeAsciiMatrix(cout, M, hdr);

  } else if (plane == "i") {

    if (idx > dims[0])
      invalidIndex(idx);
    Matrix M(dims[2]+1, dims[1]+1);
    for (int k=0; k<dims[2]; ++k)
      for (int j=0; j<dims[1]; ++j)
        M(k,j) = grid(k,j,idx);

    writeAsciiMatrix(cout, M, hdr);

  } else {
    cerr << "ERROR - unknown plane '" << plane << "'\n";
    exit(-1);
  }

}
