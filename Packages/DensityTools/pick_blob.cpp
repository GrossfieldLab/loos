/*
  pick_blob.cpp

  Given an grid-mask, a PDB, and a selection, finds the blob closest
  to ANY atom in the selection...
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
#include <algorithm>
#include <limits>

#include <DensityGrid.hpp>

using namespace std;
using namespace loos;
using namespace loos::DensityTools;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

// @cond TOOL_INTERNAL
struct Blob {
  Blob() : closest_point(0,0,0), grid_dist(numeric_limits<double>::max()),
	   real_dist(numeric_limits<double>::max()) { }

  int id;
  DensityGridpoint closest_point;
  double grid_dist;
  double real_dist;
};
// @endcond

int debug = 0;


string model_name, selection;
GCoord spot;
int picked_id;
bool use_spot = false;
bool largest = false;
double range = 0.0;

// @cond TOOLS_INTERNAL


string fullHelpMessage(void) {
  string msg =
    "SYNOPSIS\n"
    "\n"
    "\tIdentify the blob closest to a user-specified criterion\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tpick_blobs finds the blob closest to a user input.  Several methods\n"
    "of input are supported.  For instance, given a pdb and selection string\n"
    "the blob closest to the selection will be returned.  Additionally, a blob\n"
    "may be selected using its ID (see blobid).  A point within the grid may\n"
    "also be used.  A range of distances and the largest blob within the range\n"
    "are alternate criteria.\n"
    "\n"
    "The input is an integer grid (from blobid), and the output is another integer-grid\n"
    "that can be used to mask a density grid.\n"
    "\n"
    "EXAMPLES\n"
    "\tblobid --threshold 1 <foo.grid >foo_id.grid\n"
    "\tpick_blob --model foo.pdb --selection 'resid==65' < foo_id.grid > foo_picked.grid\n"
    "This example first segments the density at 1.0, and then picks the blob closest to\n"
    "any atom in residue 65 in the model.\n"
    "\n"
    "\tpick_blob --point '(13,7,3)' <foo.grid >foo_picked.grid\n"
    "This example picks the blob nearest coordinates (13,7,3) in real-space (i.e.\n"
    "Angstroms).\n"
    "\n"
    "\tpick_blob --model foo.pdb --selection 'resid==65' --range 15 <foo_id.grid >foo_picked.grid\n"
    "This example finds ALL blobs that are within 15 Angstroms of any atom in residue 65.\n"
    "\n"
    "\tpick_blob --model foo.pdb --selection 'resid==64' --range 15 --largest 1 <foo_id.grid >foo_picked.grid\n"
    "This example is as above, except that only the largest blob within 15 Angstroms is picked,\n"
    "rather than ALL blobs within 15 Angstroms.\n"
    "\n"
    "SEE ALSO\n"
    "\tblobid, gridmask\n";

  return(msg);
}


class ToolOptions : public opts::OptionsPackage {
public:
  
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("model", po::value<string>(&model_name)->default_value(""), "Select using this model (must have coords)")
      ("selection", po::value<string>(&selection)->default_value(""), "Select atoms within the PDB to find nearest blob")
      ("id", po::value<int>(&picked_id)->default_value(-1), "Select blob with this ID")
      ("point", po::value<string>(&point_spec), "Select blob closest to this point")
      ("range", po::value<double>(&range), "Select blobs that are closer than this distance")
      ("largest", po::value<bool>(&largest)->default_value(false), "Select only the largest blob that fits the distance criterion");
  }

  bool postConditions(po::variables_map& vm) {
    if (!model_name.empty()) {
      if (selection.empty()) {
        cerr << "Error: must provide a selection when using a model to select blobs\n";
        return(false);
      }
    } else if (!point_spec.empty()) {
      use_spot = true;
      istringstream iss(point_spec);
      if (!(iss >> spot)) {
        cerr << "Error: cannot parse coordinate " << point_spec << endl;
        return(false);
      }
    } else if (picked_id < 0) {
      cerr << "Error: must specify either a PDB with selection, a point, or a blob-ID to pick\n";
      return(false);
    }

    return(true);
  }

  string spot_spec;
  string point_spec;
};


void zapGrid(DensityGrid<int>& grid, const vector<int>& vals) {
  DensityGridpoint dims = grid.gridDims();

  for (int k=0; k<dims[2]; k++)
    for (int j=0; j<dims[1]; j++)
      for (int i=0; i<dims[0]; i++) {
	DensityGridpoint point(i,j,k);
	int val = grid(point);
	vector<int>::const_iterator ci = find(vals.begin(), vals.end(), val);
	if (ci == vals.end())
	  grid(point) = 0;
      }
}


int maxBlobId(const DensityGrid<int>& grid) {
  DensityGridpoint dims = grid.gridDims();
  long k = dims[0] * dims[1] * dims[2];
  int maxid = 0;

  for (long i=0; i<k; i++)
    if (grid(i) > maxid)
      maxid = grid(i);

  return(maxid);
}


vector<Blob> pickBlob(const DensityGrid<int>& grid, const vector<GCoord>& points) {
  vector<DensityGridpoint> gridded;
  vector<GCoord>::const_iterator ci;

  for (ci = points.begin(); ci != points.end(); ++ci)
    gridded.push_back(grid.gridpoint(*ci));

  int maxid = maxBlobId(grid);

  if (debug >= 1)
    cerr << boost::format("Found %d total blobs in grid.\n") % maxid;
  
  vector<Blob> blobs(maxid+1, Blob());

  DensityGridpoint dims = grid.gridDims();
  for (int k=0; k<dims[2]; k++)
    for (int j=0; j<dims[1]; j++)
      for (int i=0; i<dims[0]; i++) {
	DensityGridpoint point(i,j,k);
	int id = grid(point);
	if (!id)
	  continue;

	vector<DensityGridpoint>::iterator cj;
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
  
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  vector<GCoord> points;
  if (picked_id < 0) {
    if (use_spot)
      points.push_back(spot);
    else {
      AtomicGroup model = createSystem(model_name);
      AtomicGroup subset = selectAtoms(model, selection);
      
      for (AtomicGroup::iterator i = subset.begin(); i != subset.end(); ++i)
        points.push_back((*i)->coords());
    }
  }

  DensityGrid<int> grid;
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
