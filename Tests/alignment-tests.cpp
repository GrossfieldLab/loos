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
  
  // Suppress for easy diffs...
  //cout << invocationHeader(argc, argv) << endl;

  loos::base_generator_type& rng = loos::rng_singleton();

  // Uncommont the follong to seed the suite-wide RNG
  //loos::randomSeedRNG();


  PDB pdb(argv[1]);

  CAlphaSelector casel;
  AtomicGroup cas = pdb.select(casel);
  cout << "Found " << cas.size() << " CAlphas.\n";

  AtomicGroup casb = cas.copy();
  casb.perturbCoords(1.0);
  greal rmsd = cas.rmsd(casb);
  cout << "RMSD test = " << rmsd << endl;
  XForm W;
  W.rotate('Y', 35.0);
  casb.applyTransform(W);
  rmsd = cas.rmsd(casb);
  cout << "Rotated rmsd = " << rmsd << endl;

  GMatrix M = casb.alignOnto(cas);
  cout << "Aligned rmsd = " << cas.rmsd(casb) << endl;
  cout << M << endl;

  // -------------------------------------------------------------------------------

  cout << "*** Wide Single Alignment Sweep ***\n";
  greal x, y, z;

  for (z = -90.0; z <= 90.0; z += 30.0) {
    for (y = -90.0; y <= 90.0; y += 30.0) {
      for (x = -90.0; x <= 90.0; x += 30.0) {
	AtomicGroup casr = cas.copy();
	W.identity();
	W.rotate('x', x);
	W.rotate('y', y);
	W.rotate('z', z);
	casr.applyTransform(W);
	greal pre_rmsd = cas.rmsd(casr);
	casr.alignOnto(cas);
	greal rmsd = cas.rmsd(casr);
	cout << setw(15) << z << setw(15) << y << setw(15) << x << setw(15) << pre_rmsd << setw(15) << rmsd << endl;
      }
    }
  }


  // -------------------------------------------------------------------------------
  // Now run iteratative superpositon tests...

  boost::uniform_real<> anglemap(-45.0, 45.0);
  boost::variate_generator<loos::base_generator_type&, boost::uniform_real<> > angles(rng, anglemap);

  boost::uniform_real<> transmap(-10.0, 10.0);
  boost::variate_generator<loos::base_generator_type&, boost::uniform_real<> > translations(rng, transmap);

  vector<AtomicGroup> mols;
  for (i=0; i<maxgrpcnt; i++) {
    AtomicGroup subgroup = cas.copy();
    GCoord t(translations(), translations(), translations());
    
    subgroup.perturbCoords(2.0);

    W.identity();
    W.translate(t);
    W.rotate('z', angles());
    W.rotate('y', angles());
    W.rotate('x', angles());
	     
    subgroup.applyTransform(W);
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
