/*
  frame2pdb
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

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

// @cond TOOLS_INTERNAL


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "Extract a frame from a trajectory, writing it out as a PDB\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "Given a model, a trajectory, and a frame number, this tool will extract that\n"
    " frame and write it out as a PDB.  Optionally, a subset of the model can be \n"
    "extracted.  Any LOOS supported model and trajectory type may be used.  Note that\n"
    "frame numbers are zero-based.  Negative frame numbers are relative to the end\n"
    "of the trajectory.  Note that you will need to put '--' on the command line\n"
    "*after* any options to tell the options parse that the negative frame number\n"
    "is not another command line option.\n"
    "\n"
    "The --clear-element option is there because some build systems can produce weird\n"
    "output in the elements field of the PDB file that can cause pymol to have trouble\n"
    "rendering a protein. \n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tframe2pdb model.psf simulation.dcd 42 >frame.pdb\n"
    "Extracts the 43rd frame from the simulation.\n"
    "\n"
    "\tframe2pdb -- model.psf simulation.dcd -1 >frame.pdb\n"
    "Extracts the last frame from the simulation.\n"
    "\n"
    "\tframe2pdb --selection 'resid <= 100' model.psf simulation.dcd 13 >frame.pdb\n"
    "Extracts the 14th frame, only writing out the first 100 residues.\n";
  return(msg);
}


class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() : use_bonds(true), clear_element(false) { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("bonds", po::value<bool>(&use_bonds)->default_value(use_bonds), "Include bonds in output (if available)")
      ("clear-element", po::value<bool>(&clear_element)->default_value(clear_element), "Clear the element field in the pdb");
  }


  string print() const {
    ostringstream oss;
    oss << "use_bonds=" << use_bonds;
    return(oss.str());
  }


  bool use_bonds, clear_element;
};



// @endcond


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  opts::RequiredArguments* ropts = new opts::RequiredArguments;
  ToolOptions* topts = new ToolOptions;
  ropts->addArgument("frameno", "frame-number");

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts).add(ropts);
  if (!options.parse(argc, argv))
    exit(-1);

  if (tropts->skip)
    cerr << "WARNING- --skip is ignored by this tool\n";

  long frameno = parseStringAs<long>(ropts->value("frameno"));
  uint frame_index;
  if (frameno < 0)
    frame_index = tropts->trajectory->nframes() + frameno;
  else
    frame_index = static_cast<uint>(frameno);

  bool b = tropts->trajectory->readFrame(frame_index);
  if (!b) {
    cerr << "Could not read frame " << frame_index << " from trajectory " << tropts->traj_name << endl;
    exit(-2);
  }
  AtomicGroup subset = selectAtoms(tropts->model, sopts->selection);
  if (!topts->use_bonds)
    subset.clearBonds();

  tropts->trajectory->updateGroupCoords(subset);
  PDB pdb = PDB::fromAtomicGroup(subset);
  if (sopts->selection != "all")
    pdb.clearBonds();
  pdb.remarks().add(hdr);

  if (topts->clear_element) {
    for (auto& atom : pdb)
      atom->PDBelement(string(""));
  }
  cout << pdb;
}
