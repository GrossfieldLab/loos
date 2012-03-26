/*
  griddiff.cpp

  Subtract one grid from another (required grids to match)

*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2012, Tod D. Romo, Alan Grossfield
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

#include <DensityGrid.hpp>
#include <DensityTools.hpp>

using namespace std;
using namespace loos;
using namespace loos::DensityTools;


int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage- griddiff grid1 grid2 >grid1-grid2\n";
    exit(-1);
  }

  string hdr = invocationHeader(argc, argv);

  DensityGrid<double> grid1;
  DensityGrid<double> grid2;

  int k = 1;
  ifstream ifs1(argv[k++]);
  ifs1 >> grid1;

  ifstream ifs2(argv[k++]);
  ifs2 >> grid2;

  DensityGridpoint dims1 = grid1.gridDims();
  DensityGridpoint dims2 = grid2.gridDims();

  if (dims1 != dims2) {
    cerr << "Error- the grid sizes must match\n";
    exit(-2);
  }

  if (grid1.minCoord() != grid2.minCoord()
      || grid1.maxCoord() != grid2.maxCoord()) {
    cerr << "Error- the extents of the grids do not match\n";
    exit(-3);
  }

  for (long i = 0; i<grid1.size(); ++i)
    grid1(i) -= grid2(i);

  grid1.addMetadata(hdr);
  cout << grid1;
}
