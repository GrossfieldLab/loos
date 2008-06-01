/*
  alignment-tests.cpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Unit testing for alignment routines...

*/

#include <iostream>
#include <iomanip>
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
  unsigned int i;

  if (argc != 2) {
    cerr << "Usage- " << argv[0] << " pdbfile\n";
    exit(-1);
  }

  cout << invocationHeader(argc, argv) << endl;

  loos::base_generator_type& rng = loos::rng_singleton();

  // Uncommont the follong to seed the suite-wide RNG
  //loos::randomSeedRNG();


  PDB pdb(argv[1]);

  CAlphaSelector casel;
  AtomicGroup cas = pdb.select(casel);
  cout << "Found " << cas.size() << " CAlphas.\n";

  AtomicGroup casb = *(cas.clone());
  casb.perturbCoords(1.0);
  greal rmsd = cas.rmsd(casb);
  cout << "RMSD test = " << rmsd << endl;
  casb.xform().rotate('Y', 35.0);
  //casb.applyTransform();
  rmsd = cas.rmsd(casb);
  cout << "Rotated rmsd = " << rmsd << endl;

  GMatrix M = casb.alignOnto(cas);
  //casb.applyTransform();
  cout << "Aligned rmsd = " << cas.rmsd(casb) << endl;
  cout << M << endl;

  // -------------------------------------------------------------------------------

  cout << "*** Narrow Single Alignment Sweep ***\n";
  greal x, y, z;
  AtomicGroup casr = *(cas.clone());

  for (z = -3.0; z <= 3.0; z += 1.0) {
    for (y = -3.0; y <= 3.0; y += 1.0) {
      for (x = -3.0; x <= 3.0; x += 1.0) {
	casr.xform().identity();
	casr.xform().rotate('x', x);
	casr.xform().rotate('y', y);
	casr.xform().rotate('z', z);
	greal rad = casr.radius();
	GCoord c = casr.centroid();

	greal pre_rmsd = cas.rmsd(casr);
	casr.alignOnto(cas);
	greal rmsd = cas.rmsd(casr);
	cout << setw(15) << z << setw(15) << y << setw(15) << x << " " << c << setw(15) << rad << setw(15) << pre_rmsd << setw(15) << rmsd << endl;
      }
    }
  }

  // -------------------------------------------------------------------------------

  cout << "*** Wide Single Alignment Sweep ***\n";

  for (z = -90.0; z <= 90.0; z += 30.0) {
    for (y = -90.0; y <= 90.0; y += 30.0) {
      for (x = -90.0; x <= 90.0; x += 30.0) {
	AtomicGroup casr = *(cas.clone());
	casr.xform().identity();
	casr.xform().rotate('x', x);
	casr.xform().rotate('y', y);
	casr.xform().rotate('z', z);

	greal pre_rmsd = cas.rmsd(casr);
	casr.alignOnto(cas);
	greal rmsd = cas.rmsd(casr);
	cout << setw(15) << z << setw(15) << y << setw(15) << x << setw(15) << pre_rmsd << setw(15) << rmsd << endl;
      }
    }
  }


  // -------------------------------------------------------------------------------
  // Now run iteratative superpositon tests...

  boost::uniform_real<> uni(-45.0, 45.0);
  boost::variate_generator<loos::base_generator_type&, boost::uniform_real<> > func(rng, uni);

  vector<AtomicGroup> mols;
  for (i=0; i<maxgrpcnt; i++) {
    AtomicGroup subgroup = *(cas.clone());
    subgroup.perturbCoords(2.0);
    subgroup.xform().rotate('y', func());
    subgroup.applyTransform();
    mols.push_back(subgroup);
  }
  AtomicGroup avg = averageStructure(mols);
  cout << "Pre-aligned rmsds:\n";
  for (i=0; i<maxgrpcnt; i++)
    cout << "\t" << i << "\t" << avg.rmsd(mols[i]) << endl;


  greal final_rmsd = loos::iterativeAlignment(mols, 2.0);
  cout << "Final alignment rmsd to avg struct = " << final_rmsd << endl;

  avg = averageStructure(mols);
  for (i=0; i<maxgrpcnt; i++)
    cout << "\t" << i << "\t" << avg.rmsd(mols[i]) << endl;
}
