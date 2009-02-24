/*
  subsetter.cpp

  A general purpose tool for subsetting a trajectory.  This tool can
  be used to extract specific atoms from a trajectory or specific
  frames.  It can also be used to add periodic box information (or
  correct it) to a trajectory.  It can also be used to concatenate
  trajectories together (optionally extracting a subset of the
  concatenated trajectory).  Finally, you can center the output so
  that the centroid of the selection is the origin.  Note that the
  selection used for centering comes from the specified subset...

  The output is always in DCD format.

  Usage:
    subsetter [options] output-prefix input-model input-trajectory [input-trajectory ...]

  Examples:

  * subsetter -S10 out.dcd model.pdb traj1.dcd traj2.dcd traj3.dcd
    This concatenates the 3 trajectories together and outputs every
    10th frame to out.dcd

  * subsetter -c 'name == "CA"' out.dcd model.pdb traj1.dcd traj2.dcd traj3.dcd
    This concatenates the 3 trajectories together centering the output
    using the centroid of all c-alphas.

  * subsetter -c 'segid == "HEME"' -s '!hydrogen' out.dcd model.pdb traj.dcd
    This pulls all non-hydrogen atoms out of the trajectory and writes
    them to out.dcd, centering so that the HEME segment is at the
    origin.

  * subsetter -r 0:49 -r 150:10:300 out.dcd model.pdb traj1.dcd traj2.dcd
    This concatenates the two trajectories together, then writes out
    the first 50 frames, then frames 150 through 300 stepping by 10
    frames.  The frame indices written are of the composite trajectory.

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
#include <sstream>

using namespace std;
using namespace loos;
namespace po = boost::program_options;

// Globals...yuck...

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

bool box_override = false;
GCoord box;

string center_selection;
bool center_flag = false;





void parseOptions(int argc, char *argv[]) {

  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("verbose,v", "Verbose output")
      ("updates,u", po::value<uint>(&verbose_updates)->default_value(100), "Frequency of verbose updates")
      ("selection,s", po::value<string>(&selection)->default_value("all"), "Subset selection")
      ("stride,S", po::value<uint>(), "Step through this number of frames in each trajectory")
      ("range,r", po::value< vector<string> >(), "Frames of the DCD to use (list of Octave-style ranges)")
      ("box,b", po::value<string>(), "Override any periodic box present with this one (a,b,c)")
      ("center,c", po::value<string>(), "Recenter the trajectory using this selection (of the subset)");

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
      cerr << "Usage- " << argv[0] << " [options] output-prefix model-name trajectory-name [trajectory-name ...]\n";
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

    if (vm.count("box")) {
      string s = vm["box"].as<string>();
      istringstream is(s);
      if (!(is >> box)) {
        cerr << "ERROR - unable to convert " << s << ".  It must be in (a,b,c) format.\n";
        exit(-1);
      }
      box_override = 1;
    }

    if (vm.count("center")) {
      center_selection = vm["center"].as<string>();
      center_flag = true;
    }      

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



// Build the mapping of global frame indices into individual files,
// and into the frame number within each file...  Also returns the
// total number of frames in the composite trajectory.

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

  AtomicGroup centered;
  if (!center_selection.empty())
    centered = selectAtoms(subset, center_selection);

  uint total_frames = bindFilesToIndices(model);

  // If no frame ranges were specified, fill in all possible frames,
  // using the stride (if given)...
  if (indices.empty()) {
    stride = (!stride) ? 1 : stride;
    for (uint i=0; i<total_frames; i += stride)
      indices.push_back(i);
  }

  DCDWriter dcdout(out_name + ".dcd");
  dcdout.setTitle(hdr);

  bool first = true;  // Flag to pick off the first frame for a
                      // reference structure
  uint cnt = 0;       // Count of frames actually written

  pTraj traj;
  int current = -1;   // Track the index of the current trajectory
                      // we're reading from...

  // Iterate over all requested global-frames...
  vector<uint>::iterator vi;
  for (vi = indices.begin(); vi != indices.end(); ++vi) {

    // Have we switched to a new file??
    if (static_cast<int>(file_binding[*vi]) != current) {
      current = file_binding[*vi];
      traj = createTrajectory(traj_names[current], model);
    }

    // Read the apropriate local frame...
    traj->readFrame(local_indices[*vi]);
    traj->updateGroupCoords(model);

    // Handle centering...
    if (center_flag) {
      GCoord c = centered.centroid();
      model.translate(-c);
    }

    // Handle Periodic boundary conditions...
    if (box_override) {
      if (first && subset.isPeriodic())
        cerr << "WARNING - overriding existing periodic box.\n";
      subset.periodicBox(box);
    }

    dcdout.writeFrame(subset);

    // Pick off the first frame for the reference structure...
    if (first) {
      PDB pdb = PDB::fromAtomicGroup(subset);
      if (selection != "all")     // Strip connectivity if subsetting,
        pdb.clearBonds();         // otherwise the PDB writer will fail...
      pdb.remarks().add(hdr);
      string out_pdb_name = out_name + ".pdb";
      ofstream ofs(out_pdb_name.c_str());
      ofs << pdb;
      ofs.close();
      first = false;
    }

    ++cnt;
    if (verbose && (cnt % verbose_updates == 0))
      cerr << boost::format("Processing frame #%d (%d:%s:%d)...\n") % cnt % *vi % traj_names[current] %
        local_indices[*vi];
  }

  if (verbose)
    cerr << boost::format("Wrote %d frames to %s\n") % cnt % (out_name + ".dcd");
}
