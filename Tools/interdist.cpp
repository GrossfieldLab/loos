/*
  interdist.cpp

  (c) 2008, 2009 Tod D. Romo, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry


  Computes distances between two selections over a trajectory...

  Assumes that the trajectory has already been aligned.

  Usage:  interdist mode pdb dcd sel1 sel2 [sel3 ...]
*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008,2009 Tod Romo
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
#include <limits>
#include <boost/program_options.hpp>

using namespace std;
using namespace loos;
namespace opts = loos::OptionsFramework;




struct DistanceCalculation {
  virtual double operator()(const AtomicGroup&, const AtomicGroup&) = 0;
  virtual ~DistanceCalculation() { }
};



struct CenterDistance : public DistanceCalculation {
  double operator()(const AtomicGroup& u, const AtomicGroup& v) {
    GCoord cu = u.centroid();
    GCoord cv = v.centroid();
  
    return(cu.distance(cv));
  }
};


// Minimum distance between any member of group u vs any member of
// group v

struct MinDistance : public DistanceCalculation {
  double operator()(const AtomicGroup& u, const AtomicGroup& v) {
    double mind = numeric_limits<double>::max();

    for (AtomicGroup::const_iterator aj = v.begin(); aj != v.end(); ++aj) {
      GCoord y = (*aj)->coords();
      
      for (AtomicGroup::const_iterator ai = u.begin(); ai != u.end(); ++ai) {
        double d = y.distance2((*ai)->coords());
        if (d < mind)
          mind = d;
      }
    }
    
    return(sqrt(mind));
  }
};


// Maximum distance between any member of group u vs any member of
// group v

struct MaxDistance : public DistanceCalculation {
  double operator()(const AtomicGroup& u, const AtomicGroup& v) {
    double maxd = 0.0;

    for (AtomicGroup::const_iterator aj = v.begin(); aj != v.end(); ++aj) {
      GCoord y = (*aj)->coords();
      
      for (AtomicGroup::const_iterator ai = u.begin(); ai != u.end(); ++ai) {
        double d = y.distance2((*ai)->coords());
        if (d > maxd)
          maxd = d;
      }
    }
    
    return(sqrt(maxd));
  }
};




// @cond TOOL_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() : mode_name("center") { }

  void addGeneric(opts::po::options_description& o) {
    o.add_options()
      ("mode", opts::po::value<string>(&mode_name)->default_value(mode_name), "Calculation type (center|min|max)");
  }

  void addHidden(opts::po::options_description& o) {
    o.add_options()
      ("selection", opts::po::value< vector<string> >(&selection_names), "Selections");
  }

  void addPositional(opts::po::positional_options_description& p) {
    p.add("selection", -1);
  }

  bool check(opts::po::variables_map& map) {
    return(!selection_names.empty());
  }

  bool postConditions(opts::po::variables_map& map) {
    if (mode_name == "center")
      calc_type = new CenterDistance;
    else if (mode_name == "min")
      calc_type = new MinDistance;
    else if (mode_name == "max")
      calc_type = new MaxDistance;
    else {
      cerr << "Error- calculation mode must be 'center', 'min', or 'max'.\n";
      return(false);
    }

    return(true);
  }

  string mode_name;
  vector<string> selection_names;
  DistanceCalculation* calc_type;
};
// @endcond



int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::BasicTrajectoryOptions* tropts = new opts::BasicTrajectoryOptions;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  cout << "# " << header << endl;

  AtomicGroup model = createSystem(tropts->model_name);
  pTraj traj = createTrajectory(tropts->traj_name, model);
  vector<uint> indices = opts::assignFrameIndices(traj, tropts->frame_index_spec, tropts->skip);

  AtomicGroup src = selectAtoms(model, topts->selection_names[0]);

  cout << "# t ";
  
  vector<AtomicGroup> targets;
  for (uint i=1; i<topts->selection_names.size(); ++i) {
    AtomicGroup trg = selectAtoms(model, topts->selection_names[i]);
    targets.push_back(trg);
    cout << "d_0_" << i-1 << " ";
  }
  cout << endl;


  for (uint j=0; j<indices.size(); ++j) {
    traj->readFrame(indices[j]);
    traj->updateGroupCoords(model);

    cout << j << " ";

    for (vector<AtomicGroup>::iterator i = targets.begin(); i != targets.end(); ++i)
      cout << (*(topts->calc_type))(src, *i) << " ";
    cout << endl;
  }

}
