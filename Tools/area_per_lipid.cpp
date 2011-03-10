
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
#include <boost/program_options.hpp>

using namespace std;
using namespace loos;
namespace po = boost::program_options;


string model_name, traj_name;
string lipid_selection;
uint n_lipids;



void parseArgs(int argc, char *argv[]) {
  
  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("nlipids,n", po::value<uint>(&n_lipids)->default_value(0), "Explicit number of lipids per leaflet")
      ("headgroup,h", po::value<string>(&lipid_selection)->default_value("resname =~ 'P.GL'"), "Selection to pick lipid head groups");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model")
      ("trajectory", po::value<string>(&traj_name), "Traj");


    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("trajectory", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("trajectory"))) {
      cout << "Usage- " << argv[0] << " [options] model trajectory\n";
      cout << generic;
      exit(0);
    }

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}





int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  parseArgs(argc, argv);

  cout << "# " << hdr << endl;

  AtomicGroup model = createSystem(model_name);
  pTraj traj = createTrajectory(traj_name, model);
  if (!traj->hasPeriodicBox()) {
    cerr << "Error- trajectory has no periodicity.  Cannot compute area per lipid.\n";
    exit(-2);
  }

  // Divine how many lipids there are per leaflet...
  if (n_lipids == 0) {
    if (lipid_selection.empty()) {
      cerr << "Error- you must specify either an explicit number of lipids per leaflet or the selection to pick out the head groups\n";
      exit(-2);
    }

    traj->readFrame();
    traj->updateGroupCoords(model);
    AtomicGroup subset = selectAtoms(model, lipid_selection);
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
  cout << "# t\tArea\n";
  while (traj->readFrame()) {
    GCoord box = traj->periodicBox();
    double apl = box.x() * box.y() / n_lipids;
    cout << t++ << "\t" << apl << endl;
  }

}

