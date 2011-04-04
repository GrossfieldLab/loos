/*
  interdist.cpp

  (c) 2008, 2009, 2011 Tod D. Romo, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry

  Modified by Joshua H. Horn, Grossfield Lab
  Added functionality to consider distance only in Z-dimension

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
namespace po = boost::program_options;


bool z_only;

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

struct CenterDistanceZ : public DistanceCalculation {
  double operator()(const AtomicGroup& u, const AtomicGroup& v) {
    GCoord cu = u.centroid();
    GCoord cv = v.centroid();

    GCoord temp = cu - cv;
    return (temp.z());

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



// Globals
string model_name, traj_name;
vector<string> selection_names;
uint skip = 0;
DistanceCalculation *calc_type;


void parseOptions(int argc, char *argv[]) {
  string mode_name;

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("skip,s", po::value<uint>(&skip)->default_value(0), "Number of frames to skip at start of traj")
      ("mode,m", po::value<string>(&mode_name)->default_value("center"), "Calculation type (center|min|max|zonly)");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "model")
      ("traj", po::value<string>(&traj_name), "trajectory")
      ("selection", po::value< vector<string> >(&selection_names), "selection");
    

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);
    p.add("selection", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj") && selection_names.size() >= 2)) {
      cerr << "Usage- " << argv[0] << " [options] model trajectory sel-1 sel-2 [sel-3 ...] >output\n";
      cerr << generic;
      exit(-1);
    }

    if (mode_name == "center")
      calc_type = new CenterDistance;
    else if (mode_name == "min")
      calc_type = new MinDistance;
    else if (mode_name == "max")
      calc_type = new MaxDistance;
    else if (mode_name == "zonly")
      calc_type = new CenterDistanceZ;
    else {
      cerr << "Error- calculation mode must be either 'center', 'min', or 'max.'\n";
      exit(-1);
    }

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}



int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);
  cout << "# " << header << endl;

  AtomicGroup model = createSystem(model_name);
  pTraj traj = createTrajectory(traj_name, model);

  AtomicGroup src = selectAtoms(model, selection_names[0]);

  cout << "# t ";
  
  vector<AtomicGroup> targets;
  for (uint i=1; i<selection_names.size(); ++i) {
    AtomicGroup trg = selectAtoms(model, selection_names[i]);
    targets.push_back(trg);
    cout << "d_0_" << i-1 << " ";
  }
  cout << endl;

  uint t = skip;
  if (skip > 0)
    traj->readFrame(skip-1);

  while (traj->readFrame()) {

    traj->updateGroupCoords(model);

    cout << t++ << " ";

    vector<AtomicGroup>::iterator i;
    for (i = targets.begin(); i != targets.end(); ++i)
      cout << (*calc_type)(src, *i) << " ";
    cout << endl;
  }

}
