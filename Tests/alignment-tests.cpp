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
  vector<AtomicGroup> premols;
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
    premols.push_back(subgroup);
  }

  AtomicGroup avg = averageStructure(mols);
  if (show_results) {
    cout << "Pre-aligned rmsds:\n";
    for (i=0; i<iter_tests; i++)
      cout << "\t" << i << "\t" << avg.rmsd(mols[i]) << endl;
  }


  boost::tuple<vector<XForm>, greal, int> res = loos::iterativeAlignment(mols, 0.1);
  greal final_rmsd = boost::get<1>(res);
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

  // Now check that the composite transforms are correct...
  vector<XForm> xforms = boost::get<0>(res);
  warnings = 0;
  for (i=0; i<iter_tests; i++) {
    premols[i].applyTransform(xforms[i]);
    double irmsd = premols[i].rmsd(mols[i]);
    if (irmsd > single_rmsd_thresh) {
      ++warnings;
      cout << "WARNING - possible iterative (composite) failure at " << i << " with rmsd of " << irmsd << endl;
    }
  }
  if (warnings > 0)
    cout << "*** There were " << warnings << " possible errors detected.\n";
  else
    cout << "All composite iterative tests passed (threshold = " << single_rmsd_thresh << ")\n";

}
