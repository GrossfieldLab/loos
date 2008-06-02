/*
  pdb-tests.cpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Unit testing for PDB/AtomicGroup methods...
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


const unsigned int maxgrpcnt = MAXGRPCNT;



int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage- " << argv[0] << " pdbfile\n";
    exit(-1);
  }

  // Suppress for easy diffs...
  //cout << invocationHeader(argc, argv) << endl;

  loos::base_generator_type& rng = loos::rng_singleton();

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
  cout << "nchains = " << pdb.numberOfChains() << endl;

  AtomicGroup::BoundingBox bbox = pdb.boundingBox();
  GCoord bmin(bbox.min[0], bbox.min[1], bbox.min[2]);
  GCoord bmax(bbox.max[0], bbox.max[1], bbox.max[2]);
  cout << "Bounding box: min = " << bmin << ", max = " << bmax << endl;

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


}

