/*
  interdist.cpp

  Modified by Joshua H. Horn, Grossfield Lab
  -  Added functionality to consider distance only in Z-dimension

  Computes distances between two selections over a trajectory...

*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008 Tod Romo
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


// @cond TOOLS_INTERNAL
bool z_only;
bool segment_output = false;
double threshold = 0.0;


struct DistanceCalculation {
  virtual double operator()(const AtomicGroup&, const AtomicGroup&) = 0;
  virtual ~DistanceCalculation() { }

  void usePeriodicity(const bool flag) { _use_periodicity = flag; }
  double distance(const GCoord& u, const GCoord& v) {
  	if (_use_periodicity)
  		return(u.distance(v, _box));
  	return(u.distance(v));
  }

  void setBox(const GCoord& v) { _box = v; }

  bool _use_periodicity;
  GCoord _box;
};



struct CenterDistance : public DistanceCalculation {
  double operator()(const AtomicGroup& u, const AtomicGroup& v) {
    GCoord cu = u.centroid();
    GCoord cv = v.centroid();

    return(distance(cu, cv));

  }
};

struct CenterOfMassDistance : public DistanceCalculation {
  double operator()(const AtomicGroup& u, const AtomicGroup& v) {
    GCoord cu = u.centerOfMass();
    GCoord cv = v.centerOfMass();

    return(distance(cu, cv));

  }
};



struct CenterDistanceZ : public DistanceCalculation {
  double operator()(const AtomicGroup& u, const AtomicGroup& v) {
    GCoord cu = u.centroid();
    GCoord cv = v.centroid();

    cu.x() = cu.y() = 0.0;
    cv.x() = cv.y() = 0.0;
    return(distance(cu, cv));
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
        double d = distance(y, (*ai)->coords());
        if (d < mind)
          mind = d;
      }
    }

    return(mind);
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
        double d = distance(y, (*ai)->coords());
        if (d > maxd)
          maxd = d;
      }
    }

    return(maxd);
  }
};




class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() : mode_name("center") { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("mode", po::value<string>(&mode_name)->default_value(mode_name), "Calculation type (center|mass|min|max|zonly)")
      ("periodic", po::value<bool>(&periodic)->default_value(true), "Use periodicity in distance calculations")
      ("threshold", po::value<double>(), "Segment output using threshold distance");
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
    else if (mode_name == "zonly")
      calc_type = new CenterDistanceZ;
    else if (mode_name == "mass")
      calc_type = new CenterOfMassDistance;
    else {
      cerr << "Error- calculation mode must be either 'center', 'mass', 'min', 'max', or 'zonly'\n";
      return(false);
    }

    calc_type->usePeriodicity(periodic);

    if (map.count("threshold")) {
      threshold = map["threshold"].as<double>();
      segment_output = true;
    }

    return(true);
  }

  string help() const { return("target selection [selection ...]"); }
  string print() const {
    ostringstream oss;
    oss << boost::format("mode='%s', target='%s', selections=(%s), periodic=%d")
      % mode_name
      % target_name
      % vectorAsStringWithCommas<string>(selection_names)
      % periodic;

    if (segment_output)
      oss << boost::format("threshold=%f") % threshold;

    return(oss.str());
  }

  string mode_name;
  string target_name;
  vector<string> selection_names;
  bool periodic;
  DistanceCalculation* calc_type;
};


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "Calculate the distance between two selections over a trajectory\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "Given a model and a trajectory this tool will parse the simulation\n"
    "and return the distance between a user supplied target selection and\n"
    "any number of probe selections.\n"
    "\n"
    "There are several modes that can be selected for this tool. Each one\n"
    "specifies a different way of determining the location within the\n"
    "selection string to use in the distance calculation:\n"
    "\t center - the geometric center\n"
    "\t mass   - the center of mass\n"
    "\t min    - the minimum distance\n"
    "\t max    - the maximum distance\n"
    "\t zonly  - only the z-component\n"
    "\n"
    "\n"
    "EXAMPLE\n"
    "\n"
    "\tinterdist model.pdb traj.dcd 'name==\"CA\" && resid==133' \\\n"
    "\t    'name==\"CA\" && resid==234'\n"
    "\n"
    "Calculate the CA to CA distance between residues 133 and 234 over the\n"
    "course of trajectory traj.dcd This will print a frame number and a \n"
    "distance for each frame to stdout.\n"
    "\n"
    "\tinterdist --mode min model.pdb traj.dcd 'name==\"NE\" && resid==135'\\\n"
    "\t  'name=~\"OE.\" && resid==247'\n"
    "\n"
    "This example is similar to the first, but --mode min returns the\n"
    "minimum distance specifically.  Note the change in the second\n"
    "selection string.  Here a regular expression was supplied that\n"
    "will select either the OE1 or OE2 atom (charmm27).  The --mode min \n"
    "option will only return the distance to the closer atom.\n"
    "\n"
    "\tinterdist --mode center model.pdb traj.dcd 'segid ==\"LIG\" 'resid=15'\\\n"
    "\t    'resid=72' 'resid=13'\n"
    "\n"
    "In this example, we provide multiple selections. The resulting output\n"
    "will have 4 columns: the frame number followed by the\n"
    "centroid-to-centroid distances to residues 15, 72, and 13, in that \n"
    "order.\n"
    "\n"
    "\tinterdist --mode zonly -r 50:250  model.pdb traj.dcd 'segid==\"PROT\"' \\\n"
    "\t  'name==\"P\" && segid==\"TPE\"'\n"
    "\n"
    "Here --mode z-only indicates thatwe are only taking the z-component\n"
    "of the distance in this measurement.  the supplied range -r 50:250 \n"
    "is used to specify frames 50 to 250 for output.\n"
    "\n";
  return(msg);
    }

// @endcond


int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);


  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;

  if (topts->periodic && !traj->hasPeriodicBox()) {
    cerr << "Error- periodicity requested but trajectory is not periodic\n";
    exit(-10);
  }

  vector<uint> indices = tropts->frameList();

  AtomicGroup src = selectAtoms(model, topts->target_name);

  cout << "# " << header << endl;
  cout << "# frame ";

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

    topts->calc_type->setBox(model.periodicBox());

    for (vector<AtomicGroup>::iterator i = targets.begin(); i != targets.end(); ++i) {
      double d = (*(topts->calc_type))(src, *i);
      if (segment_output)
	cout << (d <= threshold);
      else
	cout << d;

      if (i < targets.end()-1)
	cout << '\t';
    }

    cout << endl;
  }

}
