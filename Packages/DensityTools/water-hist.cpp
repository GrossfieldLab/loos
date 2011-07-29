/*
  water-hist.cpp

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



#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include <loos.hpp>
#include <DensityGrid.hpp>
#include <water-hist-lib.hpp>
#include <DensityTools.hpp>
#include <DensityOptions.hpp>

using namespace std;
using namespace loos;
using namespace loos::DensityTools;



// @cond TOOLS_INTERAL

class WaterHistogramOptions : public opts::OptionsPackage {
public:
  WaterHistogramOptions() :
    grid_resolution(1.0),
    count_empty_voxels(false),
    rescale_density(false),
    bulk_zclip(0.0),
    bulk_zmin(0.0), bulk_zmax(0.0)
  { }

  void addGeneric(po::options_description& opts) {
    opts.add_options()
      ("gridres", po::value<double>(&grid_resolution)->default_value(grid_resolution), "Grid resolution")
      ("empty", po::value<bool>(&count_empty_voxels)->default_value(count_empty_voxels), "Count empty voxels in bulk density estimate")
      ("bulk", po::value<double>(&bulk_zclip)->default_value(bulk_zclip), "Bulk water is defined as |Z| >= k")
      ("brange", po::value<string>(), "Bulk water (--brange a,b) is defined as a <= z < b")
      ("scale", po::value<bool>(&rescale_density)->default_value(rescale_density), "Scale density by bulk estimate")
      ("clamp", po::value<string>(), "Clamp the bounding box [(x,y,z),(x,y,z)]");
  }


  bool postConditions(po::variables_map& map) {
    if (map.count("clamp")) {
      GCoord clamp_min, clamp_max;
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

    if (map.count("brange")) {
      string s = map["brange"].as<string>();
      istringstream iss(s);
      if (!(iss >> bulk_zmin)) {
        cerr << "Error- brange format is low,high\n";
        return(false);
      }
      char c;
      iss >> c;
      if (c != ',') {
        cerr << "Error- brange format is low,high\n";
        return(false);
      }
      if (!(iss >> bulk_zmax)) {
        cerr << "Error- brange format is low,high\n";
        return(false);
      }
      
    }
    return(true);
  }    

  string print() const {
    ostringstream oss;
    oss << boost::format("gridres=%f, empty=%d, bulk_zclip=%d, scale=%d, bulk_zmin=%d, bulk_zmax=%d")
      % grid_resolution
      % count_empty_voxels
      % bulk_zclip
      % rescale_density
      % bulk_zmin
      % bulk_zmax;

    if (!clamped_box.empty())
      oss << boost::format(", clamp=[%s,%s]")
        % clamped_box[0]
        % clamped_box[1];

    return(oss.str());
  }

public:
  double grid_resolution;
  bool count_empty_voxels;
  bool rescale_density;
  double bulk_zclip;
  double bulk_zmin, bulk_zmax;
  vector<GCoord> clamped_box;
};

// @endcond







int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  // Build up possible options
  opts::BasicOptions* basopts = new opts::BasicOptions;
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  opts::BasicWater* watopts = new opts::BasicWater;
  WaterHistogramOptions *xopts = new WaterHistogramOptions;

  opts::AggregateOptions options;
  options.add(basopts).add(tropts).add(watopts).add(xopts);
  if (!options.parse(argc, argv))
    exit(-1);



  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;
  vector<uint> indices = tropts->frameList();

  AtomicGroup protein = selectAtoms(model, watopts->prot_string);
  AtomicGroup water = selectAtoms(model, watopts->water_string);

  // Handle rescaling density by using a bulk-water density estimator
  BulkEstimator* est = 0;
  if (xopts->rescale_density) {
    // Double-check the clip
    traj->readFrame(indices[0]);
    traj->updateGroupCoords(protein);
    vector<GCoord> bdd = protein.boundingBox();

    if (xopts->bulk_zclip != 0.0) {
      if (xopts->bulk_zclip <= bdd[1].z())
        cerr << "***WARNING: the z-clip for bulk solvent overlaps the protein***\n";

      ZClipEstimator* myest = new ZClipEstimator(water, traj, indices, xopts->bulk_zclip, xopts->grid_resolution);
      myest->countZero(xopts->count_empty_voxels);
      est = myest;
    } else if (xopts->bulk_zmin != 0.0 || xopts->bulk_zmax != 0.0) {
      ZSliceEstimator* myest = new ZSliceEstimator(water, traj, indices, xopts->bulk_zmin, xopts->bulk_zmax, xopts->grid_resolution);
      est = myest;
    } else
      est = new NullEstimator();
  }

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
  grid.addMetadata(vectorAsStringWithCommas(options.print()));
  cout << grid;
}



