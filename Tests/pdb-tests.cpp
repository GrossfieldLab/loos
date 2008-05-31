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

  cout << invocationHeader(argc, argv) << endl;

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

  AtomicGroup casb = *(cas.clone());
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

  cas.xform().rotate('y', 45);
  cas.xform().rotate('x', 20);
  GCoord ac1 = cas[0]->coords();
  GCoord ac2 = cas.getAtomsTransformedCoord(0);

  cout << "* Transformation test:\n";
  cout << "Pre: " << ac1 << endl;
  cout << "Post: " << ac2 << endl;
  cas.applyTransform();
  GCoord ac3 = cas[0]->coords();
  cout << "Applied: " << ac3 << endl;


}

