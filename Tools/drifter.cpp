/*
  drifter.cpp

  Calculates the distance between the centroid of each frame and
  either the average (optionally of another selection) or a fixed
  point.
*/


/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo
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


using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

// @cond TOOLS_INTERNAL

struct Base {
  virtual ~Base() { }

  virtual double operator()(const GCoord& c) =0;
};



struct averageCentroid : public Base {
  averageCentroid(AtomicGroup& model, pTraj& traj) : avg(0,0,0) {
    uint n = traj->nframes();
    for (uint i=0; i<n; ++i) {
      traj->readFrame(i);
      traj->updateGroupCoords(model);
      GCoord c = model.centroid();
      avg += c;
    }
    
    avg /= n;
  }

  double operator()(const GCoord& c) {
    return(c.distance(avg));
  }

  GCoord avg;
};




struct fixedPoint : public Base {
  fixedPoint(GCoord& c) : fixed(c) { }

  double operator()(const GCoord& c) {
    return(c.distance(fixed));
  }

  GCoord fixed;
};



string fullHelpMessage(void) {
  string msg = 
    "\n"
    "SYNOPSIS\n"
    "\tCalculates the drift in the trajectory\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tdrifter calculates the distance between the centroid of a selection for each frame\n"
    "in the trajectory and either the average centroid (optionally of a different selection)\n"
    "or a fixed point.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tdrifter --selection 'segid == \"RHOD\" && name == \"CA\"' rhod.pdb rhod.dcd\n"
    "This calculates the distance between the centroid of the rhodopsin segment CA's and\n"
    "the trajectory average.\n"
    "\n"
    "\tdrifter --selection 'segid == \"TM6\" && name == \"CA\"' \\\n"
    "\t  --centroid 'name == \"CA\"' b2ar.pdb b2ar.dcd\n"
    "This calculates the distance between the centroid of TM6 (using CA's) and\n"
    "the trajectory average over all CA's\n"
    "\n"
    "\tdrifter --fixed '(4,2,2)' model.pdb model.dcd\n"
    "This calculates the distance between the centroid of the model's CA's (the default\n"
    "selection) and coordinates (4,2,2) [in Angstroms]\n";

  return(msg);
}



class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("average", "Calculate distance from selection to the average centroid (default)")
      ("centroid", po::value<string>(&centroid), "Calculate distance to the average centroid of this selection")
      ("fixed", po::value<string>(&fixed), "Calculate distance to a fixed point (x,y,z)");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("centroid='%s', fixed='%s'") % centroid % fixed;
    return(oss.str());
  }

  string centroid, fixed;
};
// @endcond






int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection* sopts = new opts::BasicSelection("name == 'CA'");
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  // Handle computation mode...
  Base* compute;

  AtomicGroup subset = selectAtoms(tropts->model, sopts->selection);

  if (! topts->centroid.empty()) {
    AtomicGroup refsub = selectAtoms(tropts->model, topts->centroid);
    compute = new averageCentroid(refsub, tropts->trajectory);
  } else if (! topts->fixed.empty()) {
    istringstream iss(topts->fixed);
    GCoord c;
    if (!(iss >> c)) {
      cerr << boost::format("Error: cannot parse %s as a corodinate.\n") % topts->fixed;
      exit(-1);
    }
    compute = new fixedPoint(c);
  } else {
    AtomicGroup subset = selectAtoms(tropts->model, sopts->selection);
    compute = new averageCentroid(subset, tropts->trajectory);
  }

  cout << "# " << hdr << endl;
  cout << "# frame d\n";
  uint t = tropts->skip;
  if (t > 0)
    tropts->trajectory->readFrame(t-1);
  else
    tropts->trajectory->rewind();
  
  while (tropts->trajectory->readFrame()) {
    tropts->trajectory->updateGroupCoords(subset);
    GCoord c = subset.centroid();
    cout << t++ << " " << (*compute)(c) << endl;
  }
}
