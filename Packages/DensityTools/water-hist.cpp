/*
  water-hist.cpp

  (c) 2009 Tod D. Romo, Grossfield Lab, URMC
*/

#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include <loos.hpp>
#include <DensityGrid.hpp>
#include <internal-water-filter.hpp>
#include <water-hist-lib.hpp>
#include <OptionsFramework.hpp>
#include <DensityOptions.hpp>

using namespace std;
using namespace loos;
using namespace loos::DensityTools;

namespace opts = loos::DensityTools::OptionsFramework;
namespace po = boost::program_options;



class WaterHistogramOptions : public opts::OptionsPackage {
public:
  WaterHistogramOptions() :
    grid_resolution(1.0),
    count_empty_voxels(false),
    rescale_density(false),
    bulk_zclip(0.0)
  { }

  void addGeneric(po::options_description& opts) {
    opts.add_options()
      ("gridres,G", po::value<double>(&grid_resolution)->default_value(grid_resolution), "Grid resolution")
      ("empty", po::value<bool>(&count_empty_voxels)->default_value(count_empty_voxels), "Count empty voxels in bulk density estimate")
      ("bulk", po::value<double>(&bulk_zclip)->default_value(bulk_zclip), "Bulk water is defined as |Z| >= k")
      ("scale", po::value<bool>(&rescale_density)->default_value(rescale_density), "Scale density by bulk estimate")
      ("clamp", po::value<string>(), "Clamp the bounding box [(x,y,z),(x,y,z)]");
  }


  bool postConditions(po::variables_map& map) {
    GCoord clamp_min, clamp_max;
    if (map.count("clamp")) {
      string s = map["clamp"].as<string>();
      stringstream ss(s);
      if (!(ss >> clamp_min)) {
        cerr << "Error- cannot parse lower bounds for box-clamp\n";
        return(false);
      }
      clamped_box.push_back(clamp_min);

      int c = ss.get();
      if (c != ',') {
        cerr << "Error- cannot parse box-clamp\n";
        return(false);
      }
      if (!(ss >> clamp_max)) {
        cerr << "Error- cannot parse upper bounds for box-clamp\n";
        return(false);
        clamped_box.push_back(clamp_max);
      }
    }
    return(true);
  }    

public:
  double grid_resolution;
  bool count_empty_voxels;
  bool rescale_density;
  double bulk_zclip;
  vector<GCoord> clamped_box;
};








int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* basopts = new opts::BasicOptions;
  opts::BasicTrajectoryOptions* trajopts = new opts::BasicTrajectoryOptions;
  opts::BasicWaterOptions* watopts = new opts::BasicWaterOptions;
  WaterHistogramOptions *xopts = new WaterHistogramOptions;

  opts::AggregateOptions options;
  options.add(basopts).add(trajopts).add(watopts).add(xopts);
  if (!options.parse(argc, argv))
    exit(-1);



  AtomicGroup model = createSystem(trajopts->model_name);
  pTraj traj = createTrajectory(trajopts->traj_name, model);
  vector<uint> indices = opts::assignFrameIndices(traj, trajopts->frame_index_spec, trajopts->skip);

  cerr << "Probe: protein = " << watopts->prot_string << endl;
  cerr << "Probe: water = " << watopts->water_string << endl;
  AtomicGroup protein = selectAtoms(model, watopts->prot_string);
  AtomicGroup water = selectAtoms(model, watopts->water_string);

  BulkEstimator* est;
  if (xopts->rescale_density) {
    // Double-check the clip
    traj->readFrame(indices[0]);
    traj->updateGroupCoords(protein);
    vector<GCoord> bdd = protein.boundingBox();
    if (xopts->bulk_zclip <= bdd[1].z())
      cerr << "***WARNING: the z-clip for bulk solvent overlaps the protein***\n";

    ZClipEstimator* myest = new ZClipEstimator(water, traj, indices, xopts->bulk_zclip, xopts->grid_resolution);
    myest->countZero(xopts->count_empty_voxels);
    est = myest;
  } else
    est = new NullEstimator();

  cerr << *est << endl;

  WaterHistogrammer wh(protein, water, est, watopts->filter_func);
  if (!xopts->clamped_box.empty()) {
    wh.setGrid(xopts->clamped_box[0]-watopts->pad, xopts->clamped_box[1]+watopts->pad, xopts->grid_resolution);
  } else
    wh.setGrid(traj, indices, xopts->grid_resolution, watopts->pad);

  wh.accumulate(traj, indices);

  long ob = wh.outOfBounds();
  if (ob)
    cerr << "***WARNING***  There were " << ob << " out of bounds waters\n";
  
  DensityGrid<double> grid = wh.grid();
  cerr << boost::format("Grid = %s x %s @ %s\n") % grid.minCoord() % grid.maxCoord() % grid.gridDims();

  if (xopts->rescale_density) {
    double d = est->bulkDensity();
    double s = est->stdDev(d);
    cerr << boost::format("Bulk density estimate = %f, std = %f\n") % d % s;
    if (xopts->rescale_density) {
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



