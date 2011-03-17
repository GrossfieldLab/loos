/*
  pick_blob.cpp


  (c) 2008 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

  Given an grid-mask, a PDB, and a selection, finds the blob closest
  to ANY atom in the selection...
*/


#include <loos.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/tuple/tuple.hpp>
#include <algorithm>
#include <sstream>
#include <limits>

#include <sgrid.hpp>

using namespace std;
using namespace loos;



struct Blob {
  Blob() : closest_point(0,0,0), grid_dist(numeric_limits<double>::max()),
	   real_dist(numeric_limits<double>::max()) { }

  int id;
  lab::SGridpoint closest_point;
  double grid_dist;
  double real_dist;
};


int debug = 0;


string pdbname, selection;
GCoord spot;
int picked_id;
bool use_spot = false;
bool largest = false;
double range = 0.0;


namespace po = boost::program_options;

void parseOptions(int argc, char *argv[]) {

  try {

    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce this help message")
      ("pdb", po::value<string>(), "Use this PDB with a selection")
      ("selection", po::value<string>(), "Select atoms within the PDB to find nearest blob")
      ("id,i", po::value<int>(&picked_id)->default_value(-1), "Select blob with this ID")
      ("point,p", po::value<string>(), "Select blob closest to this point")
      ("range,r", po::value<double>(&range), "Select blobs that are closer than this distance")
      ("largest,l", "Select only the largest blob that fits the distance criterion")
      ("verbose,v", "Verbose output")
      ("verbosity,V", po::value<int>(), "Set verbosity level");

    po::positional_options_description p;
    p.add("pdb", 1);
    p.add("selection", 2);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
	      options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
      cerr << desc;
      exit(-1);
    }

    if (vm.count("verbose"))
      debug = 1;
    if (vm.count("verbosity"))
      debug = vm["verbosity"].as<int>();

    if (vm.count("pdb")) {
      if (!vm.count("selection")) {
	cerr << "ERROR- you must give a selection along with a PDB\n";
	exit(-1);
      }

      pdbname = vm["pdb"].as<string>();
      selection = vm["selection"].as<string>();

    } else if (vm.count("point")) {

      use_spot = true;
      string spot_str = vm["point"].as<string>();
      stringstream ss(spot_str);
      if (!(ss >> spot)) {
	cerr << "ERROR- cannot parse coordinate " << spot_str << endl;
	exit(-1);
      }

    } else if (picked_id < 0) {
      cerr << "You must specify either a PDB and a selection, a point, or a blob-ID to pick\n";
      exit(-1);
    }

    if (vm.count("largest"))
      largest = true;

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
  catch(...) {
    cerr << "Exception of unknown type!\n";
  }

}


void zapGrid(lab::SGrid<int>& grid, const vector<int>& vals) {
  lab::SGridpoint dims = grid.gridDims();

  for (int k=0; k<dims[2]; k++)
    for (int j=0; j<dims[1]; j++)
      for (int i=0; i<dims[0]; i++) {
	lab::SGridpoint point(i,j,k);
	int val = grid(point);
	vector<int>::const_iterator ci = find(vals.begin(), vals.end(), val);
	if (ci == vals.end())
	  grid(point) = 0;
      }
}


int maxBlobId(const lab::SGrid<int>& grid) {
  lab::SGridpoint dims = grid.gridDims();
  long k = dims[0] * dims[1] * dims[2];
  int maxid = 0;

  for (long i=0; i<k; i++)
    if (grid(i) > maxid)
      maxid = grid(i);

  return(maxid);
}


vector<Blob> pickBlob(const lab::SGrid<int>& grid, const vector<GCoord>& points) {
  vector<lab::SGridpoint> gridded;
  vector<GCoord>::const_iterator ci;

  for (ci = points.begin(); ci != points.end(); ++ci)
    gridded.push_back(grid.gridpoint(*ci));

  int maxid = maxBlobId(grid);

  if (debug >= 1)
    cerr << boost::format("Found %d total blobs in grid.\n") % maxid;
  
  vector<Blob> blobs(maxid+1, Blob());

  lab::SGridpoint dims = grid.gridDims();
  for (int k=0; k<dims[2]; k++)
    for (int j=0; j<dims[1]; j++)
      for (int i=0; i<dims[0]; i++) {
	lab::SGridpoint point(i,j,k);
	int id = grid(point);
	if (!id)
	  continue;

	vector<lab::SGridpoint>::iterator cj;
	for (cj = gridded.begin(); cj != gridded.end(); ++cj) {
	  double d = point.distance2(*cj);
	  if (d < blobs[id].grid_dist) {
	    blobs[id].id = id;
	    blobs[id].grid_dist = d;
	    blobs[id].closest_point = point;
	    GCoord a = grid.gridToWorld(point);
	    GCoord b = grid.gridToWorld(*cj);
	    blobs[id].real_dist = a.distance2(b);
	  }
	}
      }

  if (debug > 1) {
    cerr << "* DEBUG: Blob list dump *\n";
    for (int i=0; i<=maxid; i++)
      cerr << boost::format("\tid=%d, grid_dist=%12.8g, real_dist=%12.8g\n")
	% blobs[i].id
	% blobs[i].grid_dist
	% blobs[i].real_dist;
  }

  vector<Blob> res;
  if (range == 0.0) {
    double min = numeric_limits<double>::max();
    int id = -1;
    for (int i=1; i<=maxid; i++)
      if (blobs[i].grid_dist < min) {
	min = blobs[i].grid_dist;
	id = i;
      }

    if (id > 0)
      res.push_back(blobs[id]);

  } else {
    for (int i=1; i<=maxid; i++)
      if (blobs[i].real_dist <= range)
	res.push_back(blobs[i]);
  }

  return(res);
}


int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  vector<GCoord> points;
  if (use_spot)
    points.push_back(spot);
  else {
    AtomicGroup model = createSystem(pdbname);
    AtomicGroup subset = selectAtoms(model, selection);
    
    AtomicGroup::Iterator iter(subset);
    pAtom atom;
    while (atom = iter())
      points.push_back(atom->coords());
  }

  lab::SGrid<int> grid;
  cin >> grid;

  cerr << "Read in grid with dimensions " << grid.gridDims() << endl;

  GCoord delta = grid.gridDelta();
  double voxel_volume = 1.0 / delta[0];
  for (int i=1; i<3; i++)
    voxel_volume *= (1.0 / delta[i]);
  
  if (picked_id < 0) {

    vector<Blob> picks = pickBlob(grid, points);

    if (picks.empty()) {
      cerr << "Warning - no blobs picked\n";
      exit(0);
    }

    cerr << "Picked " << picks.size() << " blobs:\n";
    vector<Blob>::const_iterator ci;
    vector<int> ids;
    for (ci = picks.begin(); ci != picks.end(); ++ci) {
      cerr << boost::format("\tid=%d, dist=%12.8g\n") % ci->id % ci->real_dist;
      ids.push_back(ci->id);
    }
      


    zapGrid(grid, ids);
    cout << grid;

  } else {
    vector<int> picks;
    picks.push_back(picked_id);
    zapGrid(grid, picks);
    cout << grid;
  }

  exit(0);
}
