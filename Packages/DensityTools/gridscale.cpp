/*
  gridscale.cpp

  Applies a constant scaling to a grid...
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

using namespace std;
using namespace loos;
using namespace loos::DensityTools;

int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  if (argc != 2) {
    cerr << "Usage- gridscale scale-value <in-grid >out-grid\n";
    cerr << "Description- scales the density values in the grid by the specified value.\n";
    cerr << "Note- the grid must be a double-precision floating point grid.\n";
    exit(0);
  }

  double konst = strtod(argv[1], 0);

  DensityGrid<double> grid;
  cin >> grid;
  grid.scale(konst);
  grid.addMetadata(hdr);
  cout << grid;
}
  
