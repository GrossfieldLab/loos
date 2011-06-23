/*
  blobid.cpp


  (c) 2008 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

  Flood-fills a grid to identify blobs...
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
class ToolOptions : public opts::OptionsPackage {

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("lower", po::value<double>(), "Sets the lower threshold for segmenting the grid")
      ("upper", po::value<double>(), "Sets the upper threshold for segmenting the grid")
      ("thresh", po::value<double>(), "Sets the threshold for segmenting the grid.");
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
  
  opts::BasicOptions* bopts = new opts::BasicOptions;
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
