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

  if (argc != 5) {
    cerr << 
      "DESCRIPTION\n\tApply a gaussian kernel convolution with a grid\n"
      "\nUSAGE\n\tgridgauss width size scaling sigma <grid >output\n"
      "Width controls the size (in grid units) of the kernel.  Size\n"
      "determines how the gaussian is mapped onto the kernel, i.e.\n"
      "-size <= x < size.  The gaussian is f(x) = exp(-0.5*(x/sigma)^2)\n"
      "and is normalized so the sum of f(x) is one, then multiplied by\n"
      "the scaling factor.\n"
      "\nEXAMPLES\n\tgridgauss 10 3 1 1 <foo.grid >foo_smoothed.grid\n"
      "This convolves the grid with a 10x10 kernel with sigma=1, and is a good\n"
      "starting point for smoothing out water density grid.\n";
    exit(0);
  }
  
  string hdr = invocationHeader(argc, argv);

  int k = 1;
  uint width = strtol(argv[k++], 0, 10);
  double scaling = strtod(argv[k++], 0);
  double normalization = strtod(argv[k++], 0);
  double sigma = strtod(argv[k++], 0);


  vector<double> kernel;
  double scaling2 = 2.0 * scaling;

  double sum = 0.0;
  for (uint i=0; i<width; ++i) {
    double x = scaling2 * i / width - scaling;
    double f = exp(-0.5*(x/sigma)*(x/sigma));
    sum += f;
    kernel.push_back(f);
  }

  double a = normalization / sum;
  for (uint i=0; i<kernel.size(); ++i)
    kernel[i] *= a;

  cerr << "Kernel (" << kernel.size() << "): ";
  copy(kernel.begin(), kernel.end(), ostream_iterator<double>(cerr, ","));
  cerr << endl;
  cerr << "normalization = " << sum << endl;

  DensityGrid<double> grid;
  cin >> grid;
  gridConvolve(grid, kernel);

  grid.addMetadata(hdr);
  cout << grid;
}


  
