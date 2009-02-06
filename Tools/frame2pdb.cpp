/*
  frame2pdb

  frame2pdb model trajectory frameno >output
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
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
#include <boost/program_options.hpp>

using namespace std;
using namespace loos;

namespace po = boost::program_options;

string model_name, traj_name, selection;
uint frameno;


void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("selection,s", po::value<string>(&selection)->default_value("all"), "Only write out this subset");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Trajectory filename")
      ("frame", po::value<uint>(&frameno), "Frame Number");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);
    p.add("frame", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj") && vm.count("frame"))) {
      cerr << "Usage- " << argv[0] << " [options] model trajectory frame# >output.pdb\n";
      cerr << generic;
      exit(-1);
    }

  }

  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  AtomicGroup model = createSystem(model_name);
  pTraj traj = createTrajectory(traj_name, model);


  bool b = traj->readFrame(frameno);
  if (!b) {
    cerr << "Could not read frame " << frameno << " from trajectory " << traj_name << endl;
    exit(-2);
  }

  AtomicGroup subset = selectAtoms(model, selection);
  traj->updateGroupCoords(subset);
  PDB pdb = PDB::fromAtomicGroup(subset);
  if (selection != "all")
    pdb.clearBonds();
  pdb.remarks().add(hdr);
  cout << pdb << endl;
}


