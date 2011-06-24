/*
  gridmask.cpp


  Given an int grid that represents picked blobs, use this as a mask
  against a double grid...
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
#include <algorithm>
#include <sstream>
#include <limits>

#include <DensityGrid.hpp>

using namespace std;
using namespace loos;
using namespace DensityTools;



int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage- gridmask <edm_grid mask_grid >masked_edm_grid\n";
    exit(-1);
  }

  DensityGrid<int> mask;
  ifstream ifs(argv[1]);
  if (!ifs) {
    cerr << "Error - cannot open " << argv[1] << " for reading.\n";
    exit(-10);
  }

  ifs >> mask;

  DensityGrid<double> data;
  cin >> data;

  DensityGridpoint dims = data.gridDims();
  DensityGridpoint ddims = mask.gridDims();
  if (dims != ddims) {
    cerr << "Error - differing dimensions between mask and density grids.\n";
    exit(-10);
  }

  long k = dims[0] * dims[1] * dims[2];
  for (long i=0; i<k; i++)
    if (!mask(i))
      data(i) = 0.0;

  cout << data;
}
