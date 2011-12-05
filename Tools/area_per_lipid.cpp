/*

  area_per_lipid.cpp

  Calculates the area-per-lipid for a trajectory

  Usage:
    area_per_lipid [options] model traj


  Requires periodic boundary information in the trajectory to
  determine leaflet area.  The selection option determines what
  residues are considered lipids (you only need to select the
  head-group).  The number of lipids can be explicitly set via the
  "--nlipids" option.
*/

/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011, Tod D. Romo
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
#include <boost/format.hpp>

using namespace std;
using namespace loos;
namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() : n_lipids(0) { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("nlipids", po::value<uint>(&n_lipids)->default_value(n_lipids), "Explicitly set the number of lipids per leaflet");
  }

  string print() const {
    ostringstream oss;
    oss << "nlipids=" << n_lipids;
    return(oss.str());
  }

  uint n_lipids;
};


// @endcond



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  
  opts::BasicOptions* basic = new opts::BasicOptions;
  opts::BasicSelection* select = new opts::BasicSelection("resname =~ 'P.GL'");
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  ToolOptions* topts = new ToolOptions;
  opts::AggregateOptions options;

  options.add(basic).add(select).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  cout << "# " << hdr << endl;

  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;
  if (!traj->hasPeriodicBox()) {
    cerr << "Error- trajectory has no periodicity.  Cannot compute area per lipid.\n";
    exit(-2);
  }

  // Divine how many lipids there are per leaflet...
  uint n_lipids = topts->n_lipids;
  if (n_lipids == 0) {
    if (select->selection.empty()) {
      cerr << "Error- you must specify either an explicit number of lipids per leaflet or the selection to pick out the head groups\n";
      exit(-2);
    }

    traj->readFrame();
    traj->updateGroupCoords(model);
    AtomicGroup subset = selectAtoms(model, select->selection);
    vector<AtomicGroup> heads = subset.splitByResidue();
    for (vector<AtomicGroup>::iterator i = heads.begin(); i != heads.end(); ++i) {
      GCoord c = (*i).centroid();
      if (c.z() > 0.0)
        ++n_lipids;
    }
    cout << "# Automatically determined " << n_lipids << " lipids per leaflet\n";
  }

  traj->rewind();
  uint t = 0;
  cout << "# frame\tArea\n";
  while (traj->readFrame()) {
    GCoord box = traj->periodicBox();
    double apl = box.x() * box.y() / n_lipids;
    cout << t++ << "\t" << apl << endl;
  }

}

