/*
  gridinfo.cpp

  Just dump the grid header info...
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

  DensityGrid<double> grid;
  if (argc == 2) {
    string fname(argv[1]);
    if (fname == "--help" || fname == "-h" || fname == "--fullhelp") {
      cerr << "Usage- gridinfo <foo.grid\n\tgridinfo foo.grid\n";
      cerr << "\nPrints out basic information about a grid\n";
      cerr << "Requires a double-precision floating point grid.\n";
      exit(-1);
    }
    ifstream ifs(argv[1]);
    if (!ifs) {
	cerr << "Error- cannot open " << argv[1] << endl;
	exit(-1);
    }
	
    ifs >> grid;
  } else
    cin >> grid;

  loos::GCoord min = grid.minCoord();
  loos::GCoord max = grid.maxCoord();
  DensityGridpoint dim = grid.gridDims();

  cout << "Grid = " << min << " x " << max << " @ " << dim << endl;
  cout << "Resolution = " << (max.x() - min.x()) / dim.x() << endl;
  cout << "Metadata: ";

  SimpleMeta meta = grid.metadata();
  if (meta.empty())
    cout << "none\n";
  else {
    cout << endl;
    copy(meta.begin(), meta.end(), ostream_iterator<SimpleMeta::value_type>(cout, "\n"));
  }
}
