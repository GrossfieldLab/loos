/*
  drifter.cpp

  Calculates the average centroid for a selection, then writes out the
  distance between the centroid of each frame's selection and the
  average
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


#include <string>
#include <sstream>
#include <vector>
#include <boost/program_options.hpp>

#include <loos.hpp>


using namespace std;
using namespace loos;
namespace po = boost::program_options;


struct Base {
  virtual ~Base() { }

  virtual double operator()(const GCoord& c, const uint t) =0;
};



struct averageCentroid : public Base {
  averageCentroid(AtomicGroup& model, pTraj traj) : avg(0,0,0) {
    uint n = traj->nframes();
    for (uint i=0; i<n; ++i) {
      traj->readFrame(i);
      traj->updateGroupCoords(model);
      GCoord c = model.centroid();
      avg += c;
    }
    
    avg /= n;
  }

  double operator()(const GCoord& c, const uint) {
    return(c.distance(avg));
  }

  GCoord avg;
};


// Note: This relies on a little trickery...  First, it assumes that
// the appropriate frame will always have been read prior to calling.
// Second, it bypass the prohibition on copying Trajectory objects by
// binding a reference to the pTraj shared pointer...  It's probably
// better not to do these sorts of things unless you understand LOOS
// reasonably well...

struct selectionCentroid : public Base {
  selectionCentroid(AtomicGroup& subset, pTraj& traj) : ref_subset(subset.copy()), traj_(traj) { }

  double operator()(const GCoord& c, const uint) {
    traj_->updateGroupCoords(ref_subset);
    GCoord ctr = ref_subset.centroid();
    return(c.distance(ctr));
  }

  AtomicGroup ref_subset;
  pTraj& traj_;
};


struct fixedPoint : public Base {
  fixedPoint(GCoord& c) : fixed(c) { }

  double operator()(const GCoord& c, const uint) {
    return(c.distance(fixed));
  }

  GCoord fixed;
};



string model_name, traj_name, selection;
Base* compute = 0;
AtomicGroup model;
pTraj traj;
AtomicGroup subset;




void parseOptions(int argc, char *argv[]) {

  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("selection,s", po::value<string>(&selection)->default_value("name == 'CA'"), "Subset selection")
      ("average,a", "Calculate distance from selection to the average centroid (default)")
      ("centroid,c", po::value<string>(), "Calculate distance to the centroid of this selection")
      ("fixed,f", po::value<string>(), "Calculate distance to a fixed point (x,y,z)");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Trajectory filenames");


    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj"))) {
      cerr << "Usage- " << argv[0] << " [options] model-name trajectory-name\n";
      cerr << generic;
      exit(-1);
    }

    model = createSystem(model_name);
    traj = createTrajectory(traj_name, model);
    subset = selectAtoms(model, selection);

    if (vm.count("centroid")) {
      string sel = vm["centroid"].as<string>();
      AtomicGroup refsub = selectAtoms(model, sel);
      compute = new selectionCentroid(refsub, traj);
    } else if (vm.count("fixed")) {
      string s = vm["fixed"].as<string>();
      stringstream ss(s);
      GCoord c;
      if (!(ss >> c)) {
        cerr << "Error - cannot parse '" << s << "' as a coordinate.  It must be (x,y,z).\n";
        exit(-1);
      }
      compute = new fixedPoint(c);
    } else
      compute = new averageCentroid(subset, traj);

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}


int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  cout << "# " << hdr << endl;
  cout << "# t d\n";
  uint n = traj->nframes();
  for (uint i=0; i<n; ++i) {
    traj->readFrame(i);
    traj->updateGroupCoords(subset);
    GCoord c = subset.centroid();
    cout << i << " " << (*compute)(c, i) << endl;
  }
}
