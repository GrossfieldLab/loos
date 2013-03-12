/*

  area_per_lipid.cpp

  Calculates the area-per-lipid for a trajectory
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
  ToolOptions() : n_lipids(0), brief(false) { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("nlipids", po::value<uint>(&n_lipids)->default_value(n_lipids), "Explicitly set the number of lipids per leaflet")
      ("brief", po::value<bool>(&brief)->default_value(brief), "Brief output (no timeseries)");
  }

  string print() const {
    ostringstream oss;
    oss << "nlipids=" << n_lipids;
    return(oss.str());
  }

  uint n_lipids;
  bool brief;
};


// @endcond


string fullHelpMessage(void)
{
string s =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "Compute the area per lipid.\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "This tool is intended for computing the area per lipid of a membrane\n"
    "system for each frame of a trajectory. To determine the number of lipids\n"
    "in a given system, a selection string must be provided, which will be\n"
    "split by residue. To override this functionality, you can specify the number\n"
    "of lipids in one leaflet with the --nlipids option.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "To calculate the area per lipid for a trajectory with PE lipid\n"
    "headgroups, you would use a command line like\n"
    "\n"
    "\tarea_per_lipid model-file traj-file --selection='resname =~\"PEGL\"'\n"
    "\n"
    "Assuming the CHARMM27-style lipid naming, the headgroup of PE lipids\n"
    "would be its own residue with the name \"PEGL\". Using selections can be\n"
    "problematic though, as it tries to guess how many are in each leaflet by\n"
    "only counting the groups with a z>0. This assumes that you have a bilayer\n"
    "and that it is centered at z=0.\n"
    "\n"
    "If you know the number of lipids in your bilayer, say 180 lipids (or 90\n"
    "in each leaflet) you can avoid the program \"guessing\" your lipid number\n"
    "and you can simplify your command line with\n"
    "\n"
    "area_per_lipid model-file traj-file --nlipids=90\n" ;

return (s);


}


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  
  opts::BasicOptions* basic = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection* select = new opts::BasicSelection("resname =~ 'P.GL'");
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  ToolOptions* topts = new ToolOptions;
  opts::AggregateOptions options;

  options.add(basic).add(select).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

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
  }

  if (tropts->skip > 0)
    traj->readFrame(tropts->skip - 1);
  else
    traj->rewind();

  vector<double> areas;
  double avg = 0.0;
  while (traj->readFrame()) {
    GCoord box = traj->periodicBox();
    double apl = box.x() * box.y() / n_lipids;
    avg += apl;
    areas.push_back(apl);
  }
  avg /= areas.size();

  double var = 0.0;
  for (uint i=0; i<areas.size(); ++i) {
    double d = areas[i] - avg;
    var += d*d;
  }
  double stdev = sqrt(var/(areas.size()-1));

  if (!topts->brief) {
    cout << "# " << hdr << endl;
    cout << "# Automatically determined " << n_lipids << " lipids per leaflet\n";
    cout << "# average = " << avg << endl;
    cout << "# stddev = " << stdev << endl;
    for (uint i=0; i<areas.size(); ++i)
      cout << i + tropts->skip << '\t' << areas[i] << endl;
  } else
    cout << n_lipids << ' ' << areas.size() << ' ' << avg << ' ' << stdev << endl;

}

