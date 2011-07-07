/*
  Gridgauss

  Apply a gaussian kernel to a grid
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
#include <boost/format.hpp>
#include <DensityGrid.hpp>
#include <GridUtils.hpp>

using namespace std;
using namespace loos;
using namespace loos::DensityTools;


int main(int argc, char *argv[]) {

  if (argc != 3) {
    cerr << "Usage- gridgauss width sigma <grid >output\n";
    exit(-1);
  }
  
  string hdr = invocationHeader(argc, argv);
  int width = atoi(argv[1]) - 1;
  double sigma = strtod(argv[2], 0);

  vector<double> kernel = gaussian1d(width, sigma);
  cerr << "Kernel (" << kernel.size() << "): ";
  copy(kernel.begin(), kernel.end(), ostream_iterator<double>(cerr, ","));
  cerr << endl;

  DensityGrid<double> grid;
  cin >> grid;
  gridConvolve(grid, kernel);

  grid.addMetadata(hdr);
  cout << grid;
}


  
