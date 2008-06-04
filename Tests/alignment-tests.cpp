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

// Configuration options...
const unsigned int iter_tests = 1000;
const double iter_perturbation = 2.0;
const double iter_rmsd_thresh = 2.0;
const double iter_final_rmsd_thresh = 1e-3;

const unsigned int single_tests = 1000;
const double single_rmsd_thresh = 1e-10;
const static bool show_results = false;



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


  boost::uniform_real<> anglemap(-180.0, 180.0);
  boost::variate_generator<loos::base_generator_type&, boost::uniform_real<> > angles(rng, anglemap);

  boost::uniform_real<> transmap(-20.0, 20.0);
  boost::variate_generator<loos::base_generator_type&, boost::uniform_real<> > translations(rng, transmap);


  // -------------------------------------------------------------------------------

  cout << "*** Single Alignment Tests ***\n";

  int warnings = 0;
  for (i=0; i<single_tests; i++) {
    AtomicGroup casr = cas.copy();
    W.identity();
    W.translate(translations(), translations(), translations());
    W.rotate('x', angles());
    W.rotate('y', angles());
    W.rotate('z', angles());
    casr.applyTransform(W);
    greal pre_rmsd = cas.rmsd(casr);
    casr.alignOnto(cas);
    double rmsd = cas.rmsd(casr);
    if (rmsd >= single_rmsd_thresh) {
      cout << "WARNING - Possible mis-alignment - pre = " << pre_rmsd << ", post = " << rmsd << endl;
      ++warnings;
    }
  }

  if (warnings > 0)
    cout << "*** There were " << warnings << " possible errors detected.\n";
  else
    cout << "All tests passed (threshold = " << single_rmsd_thresh << ")\n";



  // -------------------------------------------------------------------------------
  // Now run iteratative superpositon tests...

  vector<AtomicGroup> mols;
  for (i=0; i<iter_tests; i++) {
    AtomicGroup subgroup = cas.copy();
    GCoord t(translations(), translations(), translations());
    
    subgroup.perturbCoords(iter_perturbation);

    W.identity();
    W.translate(t);
    W.rotate('z', angles());
    W.rotate('y', angles());
    W.rotate('x', angles());
	     
    subgroup.applyTransform(W);
    mols.push_back(subgroup);
  }

  AtomicGroup avg = averageStructure(mols);
  if (show_results) {
    cout << "Pre-aligned rmsds:\n";
    for (i=0; i<iter_tests; i++)
      cout << "\t" << i << "\t" << avg.rmsd(mols[i]) << endl;
  }


  greal final_rmsd = loos::iterativeAlignment(mols, 2.0);
  if (show_results)
    cout << "Final alignment rmsd to avg struct = " << final_rmsd << endl;

  if (final_rmsd >= iter_final_rmsd_thresh)
    cout << "WARNING - final rmsd of " << final_rmsd << " is above threshold.\n";

  warnings = 0;
  avg = averageStructure(mols);
  for (i=0; i<iter_tests; i++) {
    if (show_results)
      cout << "\t" << i << "\t" << avg.rmsd(mols[i]) << endl;
    double irmsd = avg.rmsd(mols[i]);
    if (irmsd >= iter_rmsd_thresh) {
      ++warnings;
      cout << "WARNING - possible iterative failure at " << i << " with rmsd of " << irmsd << endl;
    }
  }

  if (warnings > 0)
    cout << "*** There were " << warnings << " possible errors detected.\n";
  else
    cout << "All tests passed (threshold = " << iter_final_rmsd_thresh << " : " << iter_rmsd_thresh << ")\n";
}
