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

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions;
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
  cout << pdb << endl;
}


