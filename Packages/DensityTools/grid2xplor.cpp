/*
  grid2xplor.cpp


  Converts a grid (with a number of types) into an Xplor map...
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
#include <boost/program_options.hpp>

#include <DensityGrid.hpp>
#include <xplor-edm-writer.hpp>



namespace po = boost::program_options;

using namespace std;
using namespace loos;
using namespace loos::DensityTools;



enum GridType { CHAR, INT, FLOAT, DOUBLE };

GridType gtype;
double scaling;


void parseOptions(int argc, char *argv[]) {
  string type;

  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce this help message")
      ("type", po::value<string>(&type)->default_value("double"), "Set the grid type (char, int, float, double)")
      ("scale,s", po::value<double>(&scaling)->default_value(1.0), "Scale the grid data");


    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
	      options(desc).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
      cout << desc << endl;
      exit(-1);
    }

    if (type == "double")
      gtype = DOUBLE;
    else if (type == "float")
      gtype = FLOAT;
    else if (type == "int")
      gtype = INT;
    else if (type == "char")
      gtype = CHAR;
    else {
      cerr << "Error- unknown grid type " << type << endl;
      exit(-1);
    }
  }
    catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
  catch(...) {
    cerr << "Error - exception of unknown type!\n";
  }


}



template<typename T>
DensityGrid<double> scaleGrid(DensityGrid<T>& g, const double scale) {
  DensityGridpoint dims = g.gridDims();
  long k = dims[0] * dims[1] * dims[2];
  DensityGrid<double> out(g.minCoord(), g.maxCoord(), g.gridDims());

  for (long i = 0; i<k; i++)
    out(i) = g(i) * scale;

  out.metadata(g.metadata());
  return(out);
}


int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  DensityGrid<double> edm;
  if (gtype == CHAR) {
    DensityGrid<char> grid;
    cin >> grid;
    edm = scaleGrid(grid, scaling);

  } else if (gtype == INT) {
    DensityGrid<int> grid;
    cin >> grid;
    edm = scaleGrid(grid, scaling);

  } else if (gtype == FLOAT) {
    DensityGrid<float> grid;
    cin >> grid;
    edm = scaleGrid(grid, scaling);

  } else if (gtype == DOUBLE) {
    DensityGrid<double> grid;
    cin >> grid;
    edm = scaleGrid(grid, scaling);

  } else {
    cerr << "ERROR- bad grid type internally...\n";
    exit(-2);
  }

  edm.addMetadata(header);
  GCoord min = edm.minCoord();
  GCoord max = edm.maxCoord();
  DensityGridpoint dim = edm.gridDims();
  cerr << "Read in a grid of size " << dim << endl;
  cerr << "Grid range is from " << min << " to " << max << endl;

  writeXplorEDM<double>(cout, edm);

}
