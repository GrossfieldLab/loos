/*
  gridavg

  Average grids together.  Requires that grids have the same dimensions.
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2013, Tod D. Romo, Alan Grossfield
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

  if (argc <3) {
    cerr << 
      "DESCRIPTION\n\tAverage together multiple grids\n"
      "\nUSAGE\n\tgridavg grid1 grid2 [grid3 ...] >averaged.grid\n"
      "Requires that the grids have the same dimensions.\n"
      "\nEXAMPLES\n\tgridavg water1.grid water2.grid water3.grid >water.grid\n";
    exit(0);
  }
  
  string hdr = invocationHeader(argc, argv);


  DensityGrid<double> avg;
  int k = 1;
  ifstream ifs(argv[k]);
  if (ifs.fail()) {
    cerr << "Error- cannot open " << argv[k] << endl;
    exit(-1);
  }

  ifs >> avg;
  avg.addMetadata(hdr);
  ifs.close();

  uint n = 0;
  for (++k;k < argc; ++k) {
    ifs.open(argv[k]);
    if (ifs.fail()) {
      cerr << "Error- cannot open " << argv[k] << endl;
      DensityGrid<double> grid;
      ifs >> grid;
      if (grid.gridDims() != avg.gridDims()) {
        cerr << "Error- grid in " << argv[k] << " has dimensions " << grid.gridDims() << ",\n"
             << "but was expecting it to be " << avg.gridDims() << endl;
        exit(-2);
      }
      if (grid.minCoord() != avg.minCoord() || grid.maxCoord() != avg.maxCoord())
        cerr << "Warning- real world bounds for grid in " << argv[k] << " do not match.  Proceeding anyway...\n";
      
      for (long i=0; i<grid.size(); ++i)
        avg(i) += grid(i);

      ++n;
    }
  }

  for (long i=0; i<avg.size(); ++i)
    avg(i) /= n;
  cout << avg;

}
