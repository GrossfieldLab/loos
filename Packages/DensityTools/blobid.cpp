/*
  blobid.cpp


  (c) 2008 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

  Flood-fills a grid to identify blobs...
*/


#include <loos.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/tuple/tuple.hpp>
#include <limits>
#include <ext/slist>

#include <sgrid.hpp>
#include <sgrid_utils.hpp>

namespace po = boost::program_options;

using namespace std;
using namespace loos;




boost::tuple<double, double> parseOptions(int argc, char *argv[]) {
  double low, high;

  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce this help message")
      ("lower-threshold,l", po::value<double>(), "Sets the lower threshold for segmenting the grid")
      ("upper-threshold,u", po::value<double>(), "Sets the upper threshold for segmenting the grid")
      ("threshold,t", po::value<double>(), "Sets the threshold for segmenting the grid.");

    po::positional_options_description p;
    p.add("threshold", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
	      options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
      cerr << desc << endl;
      exit(-1);
    }

    if (vm.count("threshold")) {
      low = vm["threshold"].as<double>();
      high = numeric_limits<double>::max();
    } else {
      if (!(vm.count("lower-threshold") && vm.count("upper-threshold"))) {
	cerr << "Error- you must specify either a treshold or a threshold range.\n";
	exit(-1);
      }

      low = vm["lower-threshold"].as<double>();
      high = vm["upper-threshold"].as<double>();
    }
  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
  catch(...) {
    cerr << "Error - unknown exception.\n";
    exit(-1);
  }

  boost::tuple<double, double> res(low, high);

  return(res);
}


int fill(const lab::SGridpoint seed, const int val, lab::SGrid<double>& data_grid, lab::SGrid<int>& blob_grid, const double low, const double high) {
  vector<lab::SGridpoint> blob = floodFill(seed, data_grid, val, blob_grid, lab::ThresholdRange<double>(low, high));
  return(blob.size());
}


boost::tuple<int, int, int, double> findBlobs(lab::SGrid<double>& data_grid, lab::SGrid<int>& blob_grid, const double low, const double high) {
  lab::SGridpoint dims = data_grid.gridDims();

  int id = 1;
  int min = numeric_limits<int>::max();
  int max = numeric_limits<int>::min();
  double avg = 0.0;

  for (int k=0; k<dims[2]; k++) {
    for (int j=0; j<dims[1]; j++)
      for (int i=0; i<dims[0]; i++) {
	lab::SGridpoint point(i,j,k);
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
  boost::tuple<double, double> toil = parseOptions(argc, argv);
  double lower = boost::get<0>(toil);
  double upper = boost::get<1>(toil);

  lab::SGrid<double> data;
  cin >> data;

  cerr << "Read in grid with size " << data.gridDims() << endl;

  lab::SGrid<int> blobs(data.minCoord(), data.maxCoord(), data.gridDims());
  boost::tuple<int, int, int, double> stats = findBlobs(data, blobs, lower, upper);
  cerr << boost::format("Found %d blobs in range %6.4g to %6.4g\n") % boost::get<0>(stats) % lower % upper;
  cerr << boost::format("Min blob size = %d, max blob size = %d, avg blob size = %6.4f\n")
    % boost::get<1>(stats)
    % boost::get<2>(stats)
    % boost::get<3>(stats);


  cout << blobs;
}
