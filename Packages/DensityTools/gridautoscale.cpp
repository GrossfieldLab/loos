/*
  gridautoscale.cpp

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

using namespace std;
using namespace loos;
using namespace loos::DensityTools;

typedef DensityGrid<double>   Grid;


double findPeakDensitySlice(const Grid& grid, const int nbins) {
  DensityGridpoint dims = grid.gridDims();

  int chunk_size = dims[2] / nbins;
  long volume = chunk_size * dims[1] * dims[0];

  double max_avg_density = 0.0;
  
  int kk = 0;
  for (int k = 0; k<nbins; k++) {
    // Calculate z-range...
    DensityGridpoint bottom(0,0,k*chunk_size);
    DensityGridpoint top(0,0,chunk_size*(k+1));
    
    GCoord wbottom = grid.gridToWorld(bottom);
    GCoord wtop = grid.gridToWorld(top);

    double avg = 0.0;

    for (int sk = 0; sk < chunk_size && sk+kk < dims[2]; sk++, kk++)
      for (int j=0; j<dims[1]; j++)
	for (int i=0; i<dims[0]; i++)
	  avg += grid(kk, j, i);

    avg /= volume;
    if (avg >= max_avg_density)
      max_avg_density = avg;
  }

  if (kk < dims[2]) {
    DensityGridpoint bottom(0,0,kk);
    GCoord wbottom = grid.gridToWorld(bottom);
    double avg = 0.0;

    volume = 0;
    for (; kk < dims[2]; kk++)
      for (int j=0; j<dims[1]; j++)
	for (int i=0; i<dims[0]; i++, volume++)
	  avg += grid(kk, j, i);

    DensityGridpoint top(0,0,kk);
    GCoord wtop = grid.gridToWorld(top);

    avg /= volume;
    if (avg >= max_avg_density)
      max_avg_density = avg;
  }

  return(max_avg_density);
}



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  if (argc != 1) {
    cerr << "Usage- gridautoscale <input.grid >output.grid\n";
    exit(-1);
  }


  Grid grid;
  cin >> grid;

  DensityGridpoint dims = grid.gridDims();
  double best_avg = 0.0;
  uint best_bin = 0;
  for (uint nbins = 5; nbins <= 100; ++nbins) {
    double avg = findPeakDensitySlice(grid, nbins);
    if (avg > best_avg) {
      best_avg = avg;
      best_bin = nbins;
    }
  }

  cerr << "Scaling to 1/" << best_avg << " based on " << best_bin << " bins" << endl;
  double konst = 1.0 / best_avg;
  grid.scale(konst);
  grid.addMetadata(hdr);
  ostringstream oss;
  oss << "Auto scaling = " << konst << ", bins = " << best_bin;
  grid.addMetadata(oss.str());
  cout << grid;
}
