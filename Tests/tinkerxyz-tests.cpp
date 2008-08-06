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


int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage- " << argv[0] << " TinkerXYZfile\n";
    exit(-1);
  }

  TinkerXYZ xyz(argv[1]);
  
  cout << "Read in " << xyz.size() << " atoms.\n";
  
  cout << "minId = " << xyz.minId() << endl;
  cout << "maxId = " << xyz.maxId() << endl;
  cout << "minResid = " << xyz.minResid() << endl;
  cout << "maxResid = " << xyz.maxResid() << endl;
  cout << "nresids = " << xyz.numberOfResidues() << endl;
  cout << "nsegids = " << xyz.numberOfSegids() << endl;

  vector<GCoord> bbox = xyz.boundingBox();
  cout << "Bounding box: min = " << bbox[0] << ", max = " << bbox[1] << endl;

  cout << "Centroid = " << xyz.centroid() << endl;
  cout << "Radius = " << xyz.radius() << endl;

  // -------------------------------------------------------------------------------

  vector<AtomicGroup> chains = xyz.splitByMolecule();
  cout << "Found " << chains.size() << " molecules.\n";
  unsigned int i;
  for (i = 0; i < chains.size() && i < 10; i++)
    cout << "\t" << i << "\t" << chains[i].size() << "\t" << chains[i].centroid() << endl;
  if (i < chains.size())
    cout << "...truncated...\n";

  // -------------------------------------------------------------------------------

  Parser parsed1("!(name =~ '^H')");
  KernelSelector parsed_sel1(parsed1.kernel());
  AtomicGroup grp = xyz.select(parsed_sel1);
  cout << "Found " << grp.size() << " non-hydrogen atoms via parser.\n";
  
  HeavyAtomSelector heavy;
  grp = xyz.select(heavy);
  cout << "Found " << grp.size() << " non-hydrogen atoms via HeavyAtomSelector.\n";


}

