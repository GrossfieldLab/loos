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

string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "Generate a 3D histogram of internal waters for a trajectory\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\twater-hist generates a 3-dimensional histogram for a given selection\n"
    "over the coarse of a trajectory.  This tool was originally designed\n"
    "for tracking water internal to a membrane protein, however any group\n"
    "of atoms can be substituted for \"water\" (e.g. ligand) and \"protein\".\n"
    "\n"
    "The tool first requires that you specify what atoms will be integrated.\n"
    "This is the \"water\" selection.  Next, you need to define what is considered\n"
    "\"internal\" to the protein, filtering which waters will be considered.\n"
    "This is typically done by defining a \"protein\" selection and a mode for\n"
    "filtering: axis, box, radius, or grid.  The axis mode takes the first\n"
    "principal component for the protein and picks all waters that are within\n"
    "a given radius of that axis.  The box mode uses the bounding box for the\n"
    "protein selection (i.e. any water that is within this box).  The radius\n"
    "mode picks waters that are within a given radius of any protein atom.\n"
    "Finally, the grid mode takes a grid mask and picks any waters that are\n"
    "within the masked gridpoints.\n"
    "\n"
    "The resultant density histogram can be scaled by an estimate of the\n"
    "bulk solvent density by using the --scale option along with either\n"
    "--bulk or --brange.  The former uses the average density for any Z-plane\n"
    "that is sufficiently far from 0 (i.e. |Z| >= k) whereas the latter explicitly\n"
    "takes a Z-range to average over.  Note that you must explicitly rescale\n"
    "the density by using the --scale=1 option, otherwise the estimated bulk\n"
    "solvent density will be printed only.\n"
    "\n"
    "For visualization purposes, it you are using a membrane-protein system\n"
    "and the axis mode for filtering out waters, you may end up with a plug\n"
    "of bulk water at the protein/solvent interface.  To make it more clear\n"
    "that there is a layer of bulk solvent, use the --bulked option.  This\n"
    "adds water back into the histogram based on the Z-coordinate and the\n"
    "bounding box of the protein (with an optional padding)\n"
    "\n"
    "Water-hist treats each atom as a single grid point (based on nearest)\n"
    "grid-coordinate.  This means that even though a water should cover\n"
    "multiple grid-points based on the grid resolution and water radius,\n"
    "only one grid point will be used.  For visualization then, the\n"
    "grid should be smoothed out.  This can be done via the \"gridgauss\"\n"
    "tool which convolves the grid with a gaussian kernel.  Finally,\n"
    "the grid needs to be converted to an X-Plor electron density format\n"
    "using \"grid2xplor\".  This can then be read into PyMol, VMD, or other\n"
    "visualization tools.\n"
    "\n"
    "These tools can be chained together via Unix pipes,\n"
    "   water-hist model.pdb model.dcd | gridgauss 4 2 | grid2xplor >water.xplor\n"
    "\n"
    "For more details about available options, see the help information for the\n"
    "respective tool.\n"
    "\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\twater-hist --radius=15 --bulk=25 --scale=1 b2ar.pdb b2ar.dcd | gridgauss 4 2 |\\\n"
    "\t  grid2xplor >b2ar_water.xplor\n"
    "Internal water for a GPCR with a bulk estimate, converted to Xplor EDM:\n"
    "\n"
    "\twater-hist --bulk=25 --scale=1 --bulked=20,-25:30 b2ar.pdb b2ar.dcd |\\\n"
    "\t  gridgauss 4 2 | grid2xplor >b2ar_water.xplor\n"
    "Internal water for a GPCR with the bulk solvent layer added back, converted to Xplor EDM.\n"
    "The bulk water is for any water with Z < -25 or Z > 30 and within the bounding box\n"
    "of the protein with a 20 angstrom pad.\n"
    "\n"
    "\twater-hist --radius=20 --prot='resid > 10 && resid < 25' --mode=radius |\\\n"
    "\t  gridgauss 4 2 | grid2xplor >binding.xplor\n"
    "All water within a given radius of a binding pocket, converted to Xplor EDM:\n"
    "\n"
    "\twater-hist --gridres=0.5 b2ar.pdb b2ar.dcd | gridgauss 8 4 |\\\n"
    "\t  grid2xplor >b2ar_water.xplor\n"
    "Higher resolution grid, converted to Xplor EDM:\n"
    "\n"
    "\twater-hist --prot='resname == \"PEGL\"' --water='resname === \"PEGL\"'\\\n"
    "\t  --mode=box membrane.pdb membrane.dcd >membrane.grid\n"
    "All lipid head-group density, written as LOOS grid:\n"
    "\n"
    "\twater-hist --water='resname == \"CAU\"' --mode=box b2ar.pdb b2ar.dcd >b2ar.grid\n"
    "Ligand (carazolol) Density, written as LOOS grid:\n"
    "\n"
    "NOTES\n"
    "\n"
    "When using the --bulked option, the extents of the grid are adjusted to be\n"
    "the bounding box of the protein plus the bulked pad PLUS the global pad.\n"
    "Be careful not to make the volume too large.\n"
    "\n"
    "SEE ALSO\n"
    "\tgridgauss, grid2xplor, gridstat, gridslice, blobid, pick_blob\n";

  return(msg);
}

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
      }
      clamped_box.push_back(clamp_max);
      cerr << "Warning- clamping grid to " << clamp_min << " -> " << clamp_max << endl;
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
  opts::BasicOptions* basopts = new opts::BasicOptions(fullHelpMessage());
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

  if (basopts->verbosity >= 1)
    cerr << "Filter(s): " << watopts->filter_func->name() << endl;

  // Handle rescaling density by using a bulk-water density estimator
  BulkEstimator* est;
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
  grid.addMetadata(vectorAsStringWithCommas(options.print()));
  cout << grid;
}



