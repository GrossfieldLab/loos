/*
  subsetter.cpp

  Subsets a trajectory (stripping out any atoms that don't match the given selection)
  The output is always in DCD format.

  Usage:
    subsetter [options] input-model input-trajectory output-prefix

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
#include <boost/format.hpp>
#include <boost/foreach.hpp>

using namespace std;
using namespace loos;
namespace po = boost::program_options;

vector<uint> indices;
string model_name, traj_name, prefix;
string selection;
bool verbose = false;





void parseOptions(int argc, char *argv[]) {

  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("verbose,v", "Verbose output")
      ("selection,s", po::value<string>(&selection)->default_value("all"), "Subset selection")
      ("range,r", po::value< vector<string> >(), "Frames of the DCD to use (list of Octave-style ranges)");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Trajectory filename")
      ("prefix", po::value<string>(&prefix), "Output prefix");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);
    p.add("prefix", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj") && vm.count("prefix"))) {
      cerr << "Usage- " << argv[0] << " [options] model-name trajectory-name prefix\n";
      cerr << generic;
      exit(-1);
    }

    if (vm.count("range")) {
      vector<string> ranges = vm["range"].as< vector<string> >();
      indices = parseRangeList<uint>(ranges);
    }

    if (vm.count("verbose"))
      verbose = true;

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
  AtomicGroup subset = selectAtoms(model, selection);
  subset.clearBonds();

  // Handle null-list of frames to use...
  if (indices.empty()) {
    for (uint i=0; i<traj->nframes(); ++i)
      indices.push_back(i);
  }

  // Now get ready to write the DCD...
  DCDWriter dcdout(prefix + ".dcd");
  dcdout.setTitle(hdr);

  bool first = true;
  uint cnt = 0;

  BOOST_FOREACH(uint i, indices) {
    traj->readFrame(i);
    traj->updateGroupCoords(subset);
    dcdout.writeFrame(subset);

    // Pick off the first frame for the DCD...
    if (first) {
      PDB pdb = PDB::fromAtomicGroup(subset);
      pdb.remarks().add(hdr);
      string out_pdb_name = prefix + ".pdb";
      ofstream ofs(out_pdb_name.c_str());
      ofs << pdb;
      ofs.close();
      first = false;
    }

    ++cnt;
    if (verbose && (cnt % 100 == 0))
      cerr << boost::format("Processing frame %d (%d)...\n") % cnt % i;
  }

  if (verbose)
    cerr << boost::format("Wrote %d frames to %s\n") % cnt % (prefix + ".dcd");
}
