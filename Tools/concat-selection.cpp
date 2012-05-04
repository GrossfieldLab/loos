/*
  concat-selection.cpp


  Concatenates atoms from a trajectory into a single pdb
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


#include <loos.hpp>
#include <boost/format.hpp>


using namespace std;
using namespace loos;


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tExtracts a selection from each frame of a trajectory into a PDB\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool will extract the atoms from a selection for each frame\n"
    "of a trajectory and concatenate them into one giant PDB file.  This\n"
    "can be useful for visualizing ligand locations and paths, for example.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tconcat-selection model.psf trajectory.dcd 'resname == \"CAU\"' >foo.pdb\n"
    "This extracts the residue named CAU for each frame and concatenates them\n"
    "into foo.pdb\n"
    "\n"
    "NOTES\n"
    "\tCare should be taken since the resultant PDB may be large.\n"
    ;
  return(msg);
}



int main(int argc, char *argv[]) {
  if (argc < 4) {
    cerr << "Usage: concat-selection system trajectory selection [selection...] >output.pdb\n\n"
         << fullHelpMessage();
      
    exit(-1);
  }

  string hdr = invocationHeader(argc, argv);

  AtomicGroup model = createSystem(argv[1]);
  pTraj traj = createTrajectory(argv[2], model);

  int nsegments = argc - 3;
  vector<AtomicGroup> segments(nsegments);
  vector<AtomicGroup> subsets;

  for (int i=3; i<argc; ++i) {
    AtomicGroup subset = selectAtoms(model, argv[i]);
    subsets.push_back(subset);
  }
  
  while (traj->readFrame()) {
    traj->updateGroupCoords(model);
    for (int i=0; i<nsegments; ++i) {
      for (uint j=0; j<subsets[i].size(); ++j) {
        pAtom atom(new Atom(*(subsets[i][j])));
        segments[i].append(atom);
      }
    }
  }
  
  int atomid = 1;
  AtomicGroup combined;
  for (int j=0; j<nsegments; ++j) {
    for (uint i=0; i<segments[j].size(); ++i) {
      segments[j][i]->id(atomid++);
      segments[j][i]->resid(i+1);
    }
    combined.append(segments[j]);
  }

  PDB pdb = PDB::fromAtomicGroup(combined);
  pdb.remarks().add(hdr);
  cout << pdb;
}
