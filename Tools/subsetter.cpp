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

vector<uint> indices;            // Global indices of frames to extract
vector<uint> file_binding;       // Links global frame # to the corresponding filename
vector<uint> local_indices;      // Maps global frame numbers into local frame numbers
uint stride = 0;                 // Step through frames by this amount
                                 // (unless ranges were specified)
string model_name, out_name;
vector<string> traj_names;
string selection;
bool verbose = false;
uint verbose_updates;            // Frequency of writing updates with
                                 // verbose logging..





void parseOptions(int argc, char *argv[]) {

  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("verbose,v", "Verbose output")
      ("updates,u", po::value<uint>(&verbose_updates)->default_value(100), "Frequency of verbose updates")
      ("selection,s", po::value<string>(&selection)->default_value("all"), "Subset selection")
      ("stride,S", po::value<uint>(), "Step through this number of frames in each trajectory")
      ("range,r", po::value< vector<string> >(), "Frames of the DCD to use (list of Octave-style ranges)");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value< vector<string> >(&traj_names), "Trajectory filenames")
      ("out", po::value<string>(&out_name), "Output prefix");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("out", 1);
    p.add("model", 1);
    p.add("traj", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj") && vm.count("out"))) {
      cerr << "Usage- " << argv[0] << " [options] output-prefix model-name trajectory-name\n";
      cerr << generic;
      exit(-1);
    }

    if (vm.count("range")) {
      vector<string> ranges = vm["range"].as< vector<string> >();
      indices = parseRangeList<uint>(ranges);
    }

    if (vm.count("verbose") || vm.count("updates"))
      verbose = true;

    if (vm.count("stride"))
      stride = vm["stride"].as<uint>();

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}



uint getNumberOfFrames(const string& fname, AtomicGroup& model) {
  pTraj traj = createTrajectory(fname, model);
  return(traj->nframes());
}


uint bindFilesToIndices(AtomicGroup& model) {
  uint total_frames = 0;

  for (uint j=0; j<traj_names.size(); ++j)  {
    uint n = getNumberOfFrames(traj_names[j], model);
    if (verbose)
      cout << boost::format("Trajectory \"%s\" has %d frames\n") % traj_names[j] % n;
    total_frames += n;
    for (uint i=0; i<n; ++i) {
      file_binding.push_back(j);
      local_indices.push_back(i);
    }
  }

  return(total_frames);
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  parseOptions(argc, argv);

  AtomicGroup model = createSystem(model_name);
  AtomicGroup subset = selectAtoms(model, selection);
  subset.clearBonds();

  uint total_frames = bindFilesToIndices(model);
  // Handle null-list of frames to use...
  if (indices.empty()) {
    stride = (!stride) ? 1 : stride;
    for (uint i=0; i<total_frames; i += stride)
      indices.push_back(i);
  }

  // Now get ready to write the DCD...
  DCDWriter dcdout(out_name + ".dcd");
  dcdout.setTitle(hdr);

  bool first = true;
  uint cnt = 0;

  pTraj traj;
  int current = -1;

  BOOST_FOREACH(uint i, indices) {

    // Have we switched to a new file??
    if (static_cast<int>(file_binding[i]) != current) {
      current = file_binding[i];
      traj = createTrajectory(traj_names[current], model);
    }

    traj->readFrame(local_indices[i]);
    traj->updateGroupCoords(subset);
    dcdout.writeFrame(subset);

    // Pick off the first frame for the DCD...
    if (first) {
      PDB pdb = PDB::fromAtomicGroup(subset);
      pdb.remarks().add(hdr);
      string out_pdb_name = out_name + ".pdb";
      ofstream ofs(out_pdb_name.c_str());
      ofs << pdb;
      ofs.close();
      first = false;
    }

    ++cnt;
    if (verbose && (cnt % verbose_updates == 0))
      cerr << boost::format("Processing frame #%d (%d:%s:%d)...\n") % cnt % i % traj_names[current] %
        local_indices[i];
  }

  if (verbose)
    cerr << boost::format("Wrote %d frames to %s\n") % cnt % (out_name + ".dcd");
  
}
