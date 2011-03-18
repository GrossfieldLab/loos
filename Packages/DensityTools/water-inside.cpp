/*
  water-inside.cpp


  (c) 2008 Tod D. Romo,
      Grossfield Lab,
      University of Rochester Medical and Dental School


  Applies a given set of criteria to determine whether or not a water
  is inside a protein.  A matrix is then built up where each column
  represents a timepoint in the trajectory and each row is the
  internal water state (i.e. 1 = water is inside, 0 = water is not
  inside)

  Also tracks the volume of the probe region (i.e. what's defined as
  inside, if possible) and writes out a list of atomids that describe
  which atoms go with which rows of the matrix.

*/

#include <loos.hpp>



#include <loos.hpp>
#include <cmath>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <DensityGrid.hpp>
#include <internal-water-filter.hpp>

using namespace std;
using namespace loos;
using namespace loos::DensityTools;
namespace po = boost::program_options;


typedef Math::Matrix<int, Math::ColMajor>    Matrix;


// Globals...  BOOO!
WaterFilterBase* filter_func;
double zmin, zmax;
string water_string, prot_string, model_name, traj_name, grid_name, prefix;
DensityGrid<int> the_grid;


void parseOptions(int argc, char *argv[]) {

  try {
    double pad;
    double radius;
    string mode;

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("pad,P", po::value<double>(&pad)->default_value(1.0), "Pad (for bounding box)")
      ("radius,r", po::value<double>(&radius)->default_value(10.0), "Radius (for principal axis filter)")
      ("zrange", po::value<string>(), "Clamp the volume to integrate over in Z (min:max)")
      ("water,w", po::value<string>(&water_string)->default_value("name == 'OH2'"), "Water selection")
      ("prot,p", po::value<string>(&prot_string)->default_value("name == 'CA'"), "Protein selection")
      ("grid,g", po::value<string>(), "Name of grid to use in grid-mode")
      ("mode,m", po::value<string>(&mode)->default_value("axis"), "Mode (axis|box|grid)");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Trajectory filename")
      ("prefix", po::value<string>(&prefix), "Output prefix");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);
    p.add("prefix", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj"))) {
      cerr << "Usage- water-inside [options] model-name trajectory-name prefix\n";
      cerr << generic;
      exit(-1);
    }

    // Handle modes & validation
    if (mode == "axis") {
      filter_func = new WaterFilterAxis(radius);
    } else if (mode == "box") {
      filter_func = new WaterFilterBox(pad);
    } else if (mode == "grid") {
      if (! vm.count("grid")) {
        cerr << "ERROR - you must specify a grid to use when using grid-mode\n";
        exit(-1);
      }

      string grid_name = vm["grid"].as<string>();
      ifstream ifs(grid_name.c_str());
      ifs >> the_grid;
      cerr << "Read in grid with size " << the_grid.gridDims() << endl;
      
      filter_func = new WaterFilterBlob(the_grid);

    } else {
      cerr << "ERROR - unknown mode " << mode << endl;
      exit(-1);
    }


    // Handle "decoration"
    if (vm.count("zrange")) {
      double zmin, zmax;
      string s = vm["zrange"].as<string>();
      int i = sscanf(s.c_str(), "%lf:%lf", &zmin, &zmax);
      if (i != 2) {
        cerr << boost::format("ERROR - unable to parse range '%s'\n") % s;
        exit(-1);
      }

      filter_func = new ZClippedWaterFilter(filter_func, zmin, zmax);
    }

    

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}



void writeAtomIds(const string& fname, const AtomicGroup& grp, const string& hdr) {
  ofstream ofs(fname.c_str());
  AtomicGroup::const_iterator ci;

  uint t = 0;
  ofs << "# " << hdr << endl;
  ofs << "# i\tatomid(i)\tresidue(i)\n";
  for (ci = grp.begin(); ci != grp.end(); ++ci) {
    ofs << boost::format("%u\t%d\t%s-%d\n") % t++ % (*ci)->id() % (*ci)->name() % (*ci)->resid();
  }
}


int main(int argc, char *argv[]) {
  string hdr = loos::invocationHeader(argc, argv);
  parseOptions(argc, argv);

  AtomicGroup model = loos::createSystem(model_name);
  pTraj traj = loos::createTrajectory(traj_name, model);

  AtomicGroup subset = loos::selectAtoms(model, prot_string);
  AtomicGroup waters = loos::selectAtoms(model, water_string);

  uint m = waters.size();
  uint n = traj->nframes();
  Matrix M(m,n);
  Math::Matrix<double> V(n, 1);
  cerr << boost::format("Water matrix is %d x %d.\n") % m % n;

  uint i = 0;
  cerr << "Processing - ";

  while (traj->readFrame()) {
    if (i % 250 == 0)
      cerr << ".";

    traj->updateGroupCoords(model);

    vector<int> mask = filter_func->filter(waters, subset);
    if (mask.size() != m) {
      cerr << boost::format("ERROR - returned mask has size %u but expected %u.\n") % mask.size() % m;
      exit(-10);
    }

    for (uint j=0; j<m; ++j)
      M(j,i) = mask[j];

    V(i,0) = filter_func->volume();
    ++i;
  }

  cerr << " done\n";
  writeAsciiMatrix(prefix + ".asc", M, hdr);
  writeAsciiMatrix(prefix + ".vol", V, hdr);
  writeAtomIds(prefix + ".atoms", waters, hdr);
}
