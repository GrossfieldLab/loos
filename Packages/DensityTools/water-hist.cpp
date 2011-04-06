/*
  water-hist.cpp

  (c) 2009 Tod D. Romo, Grossfield Lab, URMC
*/

#include <loos.hpp>
#include <sgrid.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include "internal_water_filter.hpp"
#include "water-lib.hpp"
#include "water-hist-lib.hpp"

using namespace std;
using namespace loos;
using namespace lab;
namespace po = boost::program_options;
namespace h2o = banal::water;

string model_name, traj_name, prot_sel, solv_sel;
double gridres, gridpad;
vector<uint> indices;
bool estimate_bulk;
double bulk_zclip;
bool rescale;
bool count_zero;

GCoord clamp_min, clamp_max;
bool clamp = false;

WaterFilter::Base* the_filter;

void fullHelp(void) {
  cerr << "Full help not available at this time...\n";
}


void parseArgs(int argc, char *argv[]) {
  string mode;
  double radius;
  double pad;

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("fullhelp", "Extended help")
      ("mode,m", po::value<string>(&mode)->default_value("axis"), "Water filter mode:\n"
       "Values:\n"
       "  axis: principal axis\n"
       "   box: protein bounding box\n")
      ("radius,r", po::value<double>(&radius)->default_value(10.0), "Radius (for principal axis filter)")
      ("boxpad,b", po::value<double>(&pad)->default_value(1.0), "Pad the bounding box for internal waters")
      ("zrange,z", po::value<string>(), "Clamp the volume to integrate over in Z (min:max)")
      ("water,w", po::value<string>(), "Add bulk water (z-slices between cutoff and bounding box) (pad,zmin:zmax)")
      ("range,R", po::value< vector<string> >(), "Frames to operate over")
      ("gridpad,P", po::value<double>(&gridpad)->default_value(0.0), "Pad the grid bounds by this much")
      ("protein,p", po::value<string>(&prot_sel)->default_value("name == 'CA'"), "Protein selection")
      ("solvent,s", po::value<string>(&solv_sel)->default_value("name == 'OH2'"), "Solvent selection")
      ("gridres,g", po::value<double>(&gridres)->default_value(0.5), "Grid resolution")
      ("estimate,e", po::value<bool>(&estimate_bulk)->default_value(false), "Estimate bulk density")
      ("countzero,c", po::value<bool>(&count_zero)->default_value(false), "Count empty voxels in the density estimate")
      ("bulkz,Z", po::value<double>(&bulk_zclip), "Only consider waters with z > this for bulk")
      ("scale,S", po::value<bool>(&rescale)->default_value(false), "Rescale density based on bulk estimator")
      ("clamp,C", po::value<string>(), "Clamp the bounding box [(x,y,z),(x,y,z)]");



    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model")
      ("traj", po::value<string>(&traj_name), "Trajectory");
    

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);
    
    if (vm.count("help") || vm.count("fulhelp") || !(vm.count("model") && vm.count("traj"))) {
      cerr << "Usage- " << argv[0] << " [options] model trajectory >output\n";
      cerr << generic;
      if (vm.count("fullhelp"))
        fullHelp();
      exit(-1);
    }
    
    if (vm.count("range")) {
      vector<string> range_spec = vm["range"].as< vector<string> >();
      indices = parseRangeList<uint>(range_spec);
    }

    
    // Handle modes & validation
    if (mode == "axis")
      the_filter = new WaterFilter::Axis(radius);
    else if (mode == "box") {
      the_filter = new WaterFilter::Box(pad);
      gridpad += pad;
    } else {
      cerr << boost::format("ERROR - unknown mode '%s'\n") % mode;
      exit(-1);
    }

    // Handle "decoration"
    if (vm.count("zrange")) {
      double zmin, zmax;
      string s = vm["zrange"].as<string>();
      int i = sscanf(s.c_str(), "%lf:%lf", &zmin, &zmax);
      if (i != 2) {
        cerr << boost::format("ERROR - unable to parse zrange '%s'\n") % s;
        exit(-1);
      }

      the_filter = new WaterFilter::ZClipped(the_filter, zmin, zmax);
    }


    if (vm.count("water")) {
      double zmin, zmax, pad;
      string s = vm["water"].as<string>();
      int i = sscanf(s.c_str(), "%lf,%lf:%lf", &pad, &zmin, &zmax);
      if (i != 3) {
        cerr << boost::format("ERROR - unable to parse bulk range '%s'\n") % s;
        exit(-1);
      }

      the_filter = new WaterFilter::Bulked(the_filter, pad, zmin, zmax);
    }


    // Verify some combinations...
    if (rescale && !estimate_bulk)
      throw(runtime_error("You must estimate bulk density in order to rescale"));

    if (vm.count("clamp")) {
      string s = vm["clamp"].as<string>();
      stringstream ss(s);
      if (!(ss >> clamp_min))
        throw(runtime_error("Cannot parse lower bounds of clamp"));
      int c = ss.get();
      if (c != ',')
        throw(runtime_error("Error parsing clamp"));
      if (!(ss >> clamp_max))
        throw(runtime_error("Cannot parse upper bounds of clamp"));
      clamp = true;
    }

    
  }
  catch(std::exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  parseArgs(argc, argv);

  AtomicGroup model = createSystem(model_name);
  pTraj traj = createTrajectory(traj_name, model);

  if (indices.empty())
    for (uint i=0; i<traj->nframes(); ++i)
      indices.push_back(i);

  AtomicGroup protein = selectAtoms(model, prot_sel);
  AtomicGroup water = selectAtoms(model, solv_sel);

  h2o::BulkEstimator* est;
  if (estimate_bulk) {
    // Double-check the clip
    traj->readFrame(indices[0]);
    traj->updateGroupCoords(protein);
    vector<GCoord> bdd = protein.boundingBox();
    if (bulk_zclip <= bdd[1].z())
      cerr << "***WARNING: the z-clip for bulk solvent overlaps the protein***\n";

    h2o::ZClipEstimator* myest = new h2o::ZClipEstimator(water, traj, indices, bulk_zclip, gridres);
    myest->countZero(count_zero);
    est = myest;
  } else
    est = new h2o::NullEstimator();

  cerr << *est << endl;

  h2o::WaterHistogrammer wh(protein, water, est, the_filter);
  if (clamp) {
    clamp_min -= gridpad;
    clamp_max += gridpad;
    wh.setGrid(clamp_min, clamp_max, gridres);
  } else
    wh.setGrid(traj, indices, gridres, gridpad);

  wh.accumulate(traj, indices);

  long ob = wh.outOfBounds();
  if (ob)
    cerr << "***WARNING***  There were " << ob << " out of bounds waters\n";
  
  SGrid<double> grid = wh.grid();
  cerr << boost::format("Grid = %s x %s @ %s\n") % grid.minCoord() % grid.maxCoord() % grid.gridDims();

  if (estimate_bulk) {
    double d = est->bulkDensity();
    double s = est->stdDev(d);
    cerr << boost::format("Bulk density estimate = %f, std = %f\n") % d % s;
    if (rescale) {
      cerr << "Rescaling grid by bulk estimate...\n";
      grid.scale(1.0 / d);
      stringstream ss;
      ss << boost::format("water-hist: bulk density estimate = %f, std = %f") % d % s;
      grid.addMetadata(ss.str());
    }
  }

  grid.addMetadata(hdr);
  cout << grid;
}



