/*
  convert2pdb


  Converts a LOOS-supported format to a PDB (so long as coordinates
  are present)

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
    "\tConvert any LOOS model file to a PDB\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tReads in any LOOS model file and writes it to stdout as a PDB.  A subset\n"
    "of the model may be selected.  As not all formats contain coordinates,\n"
    "these may be taken from another source by using the --coordinates option.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tconvert2pdb model.gro >model.pdb\n"
    "Converts a GROMACS .gro file to a PDB\n"
    "\n"
    "\tconvert2pdb --coordinates model.rst model.prmtop >model.pdb\n"
    "Converts an AMBER PRMTOP file (taking coordinates from the RST file).\n"
    "\n"
    "\tconvert2pdb --selection 'name == \"CA\"' model.gro >model.pdb\n"
    "Converts a GROMACS .gro file to a PDB, only writing out the alpha-carbons.\n"
    "\n"
    ;

  return(msg);
}

// @endcond


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions *bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::ModelWithCoords* mwcopts = new opts::ModelWithCoords;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mwcopts);
  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup subset = selectAtoms(mwcopts->model, sopts->selection);
  PDB pdb = PDB::fromAtomicGroup(subset);
  pdb.remarks().add(hdr);
  cout << pdb;
}
