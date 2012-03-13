/*
  grid2ascii

  Convert a grid into a serialized ascii representation
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo, Alan Grossfield
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

#include <boost/format.hpp>

using namespace std;
using namespace loos;
using namespace loos::DensityTools;

int main(int argc, char *argv[]) {
  DensityGrid<double> grid;

  if (argc != 1) {
    cerr << "Usage- grid2ascii <foo.grid >foo.asc\n"
      "\n"
      "Converts a LOOS grid to an ASCII representation.  Requires a double precision\n"
      "floating point grid.\n";
    exit(-1);
  }
  
  cin >> grid;
  DensityGridpoint dim = grid.gridDims();
  GCoord min = grid.minCoord();
  GCoord max = grid.maxCoord();

  cout << boost::format("Read in grid of size %s\n") % dim;
  cout << boost::format("Grid range from %s x %s\n") % min % max;

  for (int k=0; k<dim.z(); ++k)
    for (int j=0; j<dim.y(); ++j)
      for (int i=0; i<dim.x(); ++i)
        cout << boost::format("(%d,%d,%d) = %f\n") % k % j % i % grid(k,j,i);
}
