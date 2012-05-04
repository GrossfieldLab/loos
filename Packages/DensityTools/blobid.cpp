/*
  blobid.cpp

  Flood-fills a grid to identify blobs...
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
#include <limits>
#include <ext/slist>

#include <DensityTools.hpp>
#include <DensityGrid.hpp>
#include <GridUtils.hpp>

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

using namespace std;
using namespace loos;
using namespace loos::DensityTools;

double lower, upper;

// @cond TOOLS_INTERNAL


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "\tIdentify blobs in a density grid.\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tblobid identifies blobs by density values either in a range or above a threshold.\n"
    "An edm grid (see for example water-hist) is expected for input.\n"
    "Blobid then uses a flood-fill to determine how many separate blobs\n"
    "meet the threshold/range criteria.  A new grid is then written out\n"
    "which identifies the separate blobs.\n"
    "\nEXAMPLES\n"
    "\tblobid --threshold 1 <foo.grid >foo_id.grid\n"
    "Here we include all blobs above the threshold 1.  foo_grid is a density\n"
    "grid that has been created previously.  For example a smoothed water \n"
    "histogram grid may be used: \n"
    "\twater-hist --radius=15 --bulk=25 --scale=1 b2ar.pdb b2ar.dcd |\\\n"
    "\t  grid2gauss 4 2 > foo_grid\n"
    "The resulting blobs are then written to the grid \"foo_id\"\n"
    "\n\n";

  return(msg);
}

class ToolOptions : public opts::OptionsPackage {

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("lower", po::value<double>(), "Sets the lower threshold for segmenting the grid")
      ("upper", po::value<double>(), "Sets the upper threshold for segmenting the grid")
      ("threshold", po::value<double>(), "Sets the threshold for segmenting the grid.");
  }

  bool postConditions(po::variables_map& vm) {
    if (vm.count("threshold")) {
      lower = vm["threshold"].as<double>();
      upper = numeric_limits<double>::max();
    } else {
      if (!(vm.count("lower-threshold") && vm.count("upper-threshold"))) {
	cerr << "Error- you must specify either a treshold or a threshold range.\n";
        return(false);
      }

      lower = vm["lower-threshold"].as<double>();
      upper = vm["upper-threshold"].as<double>();
    }
    return(true);
  }

  string print() const {
    ostringstream oss;

    oss << boost::format("lower=%f, upper=%f") % lower % upper;
    return(oss.str());
  }

};
// @endcond




int fill(const DensityGridpoint seed, const int val, DensityGrid<double>& data_grid, DensityGrid<int>& blob_grid, const double low, const double high) {
  vector<DensityGridpoint> blob = floodFill(seed, data_grid, val, blob_grid, ThresholdRange<double>(low, high));
  return(blob.size());
}


boost::tuple<int, int, int, double> findBlobs(DensityGrid<double>& data_grid, DensityGrid<int>& blob_grid, const double low, const double high) {
  DensityGridpoint dims = data_grid.gridDims();

  int id = 1;
  int min = numeric_limits<int>::max();
  int max = numeric_limits<int>::min();
  double avg = 0.0;

  for (int k=0; k<dims[2]; k++) {
    for (int j=0; j<dims[1]; j++)
      for (int i=0; i<dims[0]; i++) {
	DensityGridpoint point(i,j,k);
	if (blob_grid(point) == 0 && (data_grid(point) >= low && data_grid(point) <= high)) {
	  int n = fill(point, id++, data_grid, blob_grid, low, high);
	  if (n < min)
	    min = n;
	  if (n > max)
	    max = n;
	  avg += n;
	}
      }
  }

  avg /= (id-1);
  boost::tuple<int, int, int, double> res(id-1, min, max, avg);
  return(res);
}


int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);
  
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);



  DensityGrid<double> data;
  cin >> data;

  cerr << "Read in grid with size " << data.gridDims() << endl;

  DensityGrid<int> blobs(data.minCoord(), data.maxCoord(), data.gridDims());
  boost::tuple<int, int, int, double> stats = findBlobs(data, blobs, lower, upper);
  cerr << boost::format("Found %d blobs in range %6.4g to %6.4g\n") % boost::get<0>(stats) % lower % upper;
  cerr << boost::format("Min blob size = %d, max blob size = %d, avg blob size = %6.4f\n")
    % boost::get<1>(stats)
    % boost::get<2>(stats)
    % boost::get<3>(stats);


  cout << blobs;
}
