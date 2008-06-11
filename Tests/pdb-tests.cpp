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


#include <boost/random.hpp>
#include <ctime>


#include <loos.hpp>
#include <pdb.hpp>
#include <Kernel.hpp>
#include <Parser.hpp>
#include <Selectors.hpp>
#include <ensembles.hpp>


#if !defined(MAXGRPCNT)
#define MAXGRPCNT 10
#endif


const static unsigned int maxgrpcnt = MAXGRPCNT;
static loos::base_generator_type& rng = loos::rng_singleton();

void test_findById(PDB& pdb) {
  int i;

  cout << "Testing findById()...";
  int maxid = pdb.maxId();
  int minid = pdb.minId();
  boost::uniform_int<> atomid_map(-maxid, maxid);
  boost::variate_generator<loos::base_generator_type&, boost::uniform_int<> > random_ids(rng, atomid_map);
  for (i = 0; i < 25000; i++) {
    int id = random_ids();
    pAtom pa = pdb.findById(id);
    if (pa == 0) {
      if (id >= minid && id <= maxid) {
	cout << "***ERROR***  Couldn't find an atomid that we expected to find.\n";
	exit(-10);
      }
    }
    else
      if (pa->id() != id) {
	cout << "***ERROR*** Atom found was not the atom expected...\n";
	exit(-10);
      }
    
  }

  cout << "passed\n";
}


void test_selections(PDB& pdb) {
  
  cout << "Testing selections...\n";
  cout << "Total input size = " << pdb.size() << endl;
  CAlphaSelector casel;
  AtomicGroup grp = pdb.select(casel);
  cout << "CASelector = " << grp.size() << " @ " << grp.centroid() << endl;

  BackboneSelector bbsel;
  grp = pdb.select(bbsel);
  cout << "BackboneSelector = " << grp.size() << " @ " << grp.centroid() << endl;

  SegidSelector segsel("BULK");
  grp = pdb.select(segsel);
  cout << "SegidSelector(BULK) = " << grp.size() << " @ " << grp.centroid() << endl;

  ResidRangeSelector rrange(10,20);
  grp = pdb.select(rrange);
  cout << "ResidRangeSelector(10,20) = " << grp.size() << " @ " << grp.centroid() << endl;

  HydrogenSelector hsel;
  grp = pdb.select(hsel);
  cout << "HyrdogenSelector = " << grp.size() << " @ " << grp.centroid() << endl;

  HeavyAtomSelector hesel;
  grp = pdb.select(hesel);
  cout << "HeavyAtomSelector = " << grp.size() << " @ " << grp.centroid() << endl;

  SolventSelector solsel;
  grp = pdb.select(solsel);
  cout << "SolventSelector = " << grp.size() << " @ " << grp.centroid() << endl;

  NotSelector notsolsel(solsel);
  grp = pdb.select(notsolsel);
  cout << "NotSelector(SolventSelector) = " << grp.size() << " @ " << grp.centroid() << endl;
}


int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage- " << argv[0] << " pdbfile\n";
    exit(-1);
  }

  // Suppress for easy diffs...
  //cout << invocationHeader(argc, argv) << endl;

  //loos::base_generator_type& rng = loos::rng_singleton();

  // Uncommont the follong to seed the suite-wide RNG
  //loos::randomSeedRNG();


  PDB pdb(argv[1]);
  
  cout << "Read in " << pdb.size() << " atoms.\n";
  if (pdb.isPeriodic())
    cout << "Periodic boundary conditions: " << pdb.periodicBox() << endl;
  
  cout << "minId = " << pdb.minId() << endl;
  cout << "maxId = " << pdb.maxId() << endl;
  cout << "minResid = " << pdb.minResid() << endl;
  cout << "maxResid = " << pdb.maxResid() << endl;
  cout << "nresids = " << pdb.numberOfResidues() << endl;
  cout << "nsegids = " << pdb.numberOfSegids() << endl;

  vector<GCoord> bbox = pdb.boundingBox();
  cout << "Bounding box: min = " << bbox[0] << ", max = " << bbox[1] << endl;

  cout << "Centroid = " << pdb.centroid() << endl;
  cout << "Radius = " << pdb.radius() << endl;

  // -------------------------------------------------------------------------------

  CAlphaSelector casel;
  AtomicGroup cas = pdb.select(casel);
  cout << "Found " << cas.size() << " CAs\n";

  // -------------------------------------------------------------------------------

  AtomicGroup casb = cas.copy();
  casb.perturbCoords(5.0);
  greal rmsd = cas.rmsd(casb);
  cout << "RMSD test = " << rmsd << endl;

  // -------------------------------------------------------------------------------

  vector<AtomicGroup> chains = pdb.splitByUniqueSegid();
  cout << "Found " << chains.size() << " unique segids.\n";
  unsigned int i;
  for (i = 0; i < chains.size() && i < 10; i++)
    cout << "\t" << i << "\t" << chains[i].size() << "\t" << chains[i].centroid() << endl;
  if (i < chains.size())
    cout << "...truncated...\n";

  AtomicGroup grp = cas.subset(0, 3);
  cout << "* First 3 cas *\n" << grp << endl;

  grp = cas.subset(-3, 3);
  cout << "* Last 3 cas *\n" << grp << endl;

  // -------------------------------------------------------------------------------

  Parser parsed1("!(name =~ '^H')");
  KernelSelector parsed_sel1(parsed1.kernel());
  grp = pdb.select(parsed_sel1);
  cout << "Found " << grp.size() << " non-hydrogen atoms via parser.\n";
  
  HeavyAtomSelector heavy;
  grp = pdb.select(heavy);
  cout << "Found " << grp.size() << " non-hydrogen atoms via HeavyAtomSelector.\n";

  cout << "Residue for third CA:\n";
  grp = pdb.getResidue(cas[2]);
  cout << grp << endl;

  // -------------------------------------------------------------------------------

  XForm W;
  W.rotate('y', 45);
  W.rotate('x', 20);
  GCoord ac1 = cas[0]->coords();
  GCoord ac2 = W.transform(ac1);

  cout << "* Transformation test:\n";
  cout << "Pre: " << ac1 << endl;
  cout << "Post: " << ac2 << endl;

  // -------------------------------------------------------------------------------
  test_findById(pdb);
  test_selections(pdb);


}

