/*
  pdbselect.cpp

  Takes a PDB and a selection string.  Parses the selection, then
  applies it to the PDB and writes the output to stdout.  This tool is
  used maily for checking your selection strings to make sure you're
  actually selecting what you intend to select...
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



#include <pdb.hpp>
#include <Parser.hpp>
#include <Selectors.hpp>


int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);

  if (argc != 3) {
    cerr << "Usage- " << argv[0] << " <selection string> <pdb file>\n";
    exit(-1);
  }

  PDB pdb(argv[2]);
  Parser parsed(argv[1]);
  KernelSelector parsed_selector(parsed.kernel());
  
  AtomicGroup subset = pdb.select(parsed_selector);

  cerr << "You selected " << subset.size() << " atoms out of " << pdb.size() << endl;
  
  PDB output = PDB::fromAtomicGroup(subset);
  output.remarks().add(header);
  cout << output;
}
