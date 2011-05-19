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

using namespace std;
using namespace loos;
namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;




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

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("mode", po::value<string>(&mode_name)->default_value(mode_name), "Calculation type (center|min|max)");
  }

  void addHidden(po::options_description& o) {
    o.add_options()
      ("target", po::value<string>(&target_name), "Target")
      ("selection", po::value< vector<string> >(&selection_names), "Selections");
  }

  void addPositional(po::positional_options_description& p) {
    p.add("target", 1);
    p.add("selection", -1);
  }

  bool check(po::variables_map& map) {
    return(selection_names.empty() || target_name.empty());
  }

  bool postConditions(po::variables_map& map) {
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

  string help() const { return("target selection [selection ...]"); }
  string print() const {
    ostringstream oss;
    oss << boost::format("mode='%s', target='%s', selections=(%s)") 
      % mode_name
      % target_name
      % vectorAsStringWithCommas<string>(selection_names);

    return(oss.str());
  }

  string mode_name;
  string target_name;
  vector<string> selection_names;
  DistanceCalculation* calc_type;
};
// @endcond



int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  cout << "# " << header << endl;

  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;
  vector<uint> indices = tropts->frameList();

  AtomicGroup src = selectAtoms(model, topts->target_name);

  cout << "# t ";
  
  vector<AtomicGroup> targets;
  for (uint i=0; i<topts->selection_names.size(); ++i) {
    AtomicGroup trg = selectAtoms(model, topts->selection_names[i]);
    targets.push_back(trg);
    cout << "d_0_" << i << " ";
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
