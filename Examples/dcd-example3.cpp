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


// Select for non-solvent atoms...
struct NotSolvSelector : public AtomSelector {
  bool operator()(const pAtom& atom) const {
    return(!(atom->segid() == "SOLV" || atom->segid() == "BULK"));
  }
};



int main(int argc, char *argv[]) {

  // First, read in the PSF...
  PSF psf(argv[1]);

  // Extract the non-solvent atoms...
  NotSolvSelector ns;
  AtomicGroup nonsolv = psf.select(ns);
  cout << "Found " << nonsolv.size() << " non-solvent atoms.\n";

  DCD dcd(argv[2]);

  int frameno = 0;
  while (dcd.readFrame()) {
    dcd.updateGroupCoords(nonsolv);
    cout << setw(6) << frameno++ << " = " << nonsolv.centroid() << endl;
  }

}
