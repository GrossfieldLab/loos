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
    "frame numbers are zero-based.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tframe2pdb model.psf simulation.dcd 42 >frame.pdb\n"
    "Extracts the 43rd frame from the simulation.\n"
    "\n"
    "\tframe2pdb --selection 'resid <= 100' model.psf simulation.dcd 13 >frame.pdb\n"
    "Extracts the 14th frame, only writing out the first 100 residues.\n";
  return(msg);
}


// @endcond


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  opts::RequiredArguments* ropts = new opts::RequiredArguments;
  ropts->addArgument("frameno", "frame-number");
  
  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(ropts);
  if (!options.parse(argc, argv))
    exit(-1);

  if (tropts->skip)
    cerr << "WARNING- --skip is ignored by this tool\n";

  uint frameno = parseStringAs<uint>(ropts->value("frameno"));
  bool b = tropts->trajectory->readFrame(frameno);
  if (!b) {
    cerr << "Could not read frame " << frameno << " from trajectory " << tropts->traj_name << endl;
    exit(-2);
  }
  AtomicGroup subset = selectAtoms(tropts->model, sopts->selection);

  tropts->trajectory->updateGroupCoords(subset);
  PDB pdb = PDB::fromAtomicGroup(subset);
  if (sopts->selection != "all")
    pdb.clearBonds();
  pdb.remarks().add(hdr);
  cout << pdb;
}


