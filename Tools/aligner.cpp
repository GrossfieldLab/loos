/*
  aligner.cpp


  Aligns structures in a trajectory...

  Usage:

    aligner [options] structure-file trajectory-file output-prefix

  Notes:

  Takes two selections.  The first is the subset of atoms that will
  be used for the alignment.  The second is the subset of atoms that
  will then be transformed by the alignment and written out.  The
  alignment scheme is to calculate an average structure, align all
  frames against this average, then compute a new average.  This is
  iterated until the difference in average structures is below the
  specified tolerance.

  The output will ALWAYS be a DCD!

  Aligner will cache the entire alignment selection in memory, so
  beware potential memory issues...

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


#include <loos.hpp>

#include <boost/program_options.hpp>


using namespace std;
using namespace loos;

namespace po = boost::program_options;

string model_name, traj_name, prefix;
string alignment_string;
string transform_string;

string reference_name;
string reference_sel;

greal alignment_tol;

int maxiter;
bool center;



void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help,h", "Produce this help message")
      ("align,a", po::value<string>(&alignment_string)->default_value("name == 'CA'"),"Align using this selection")
      ("transform,t", po::value<string>(&transform_string)->default_value("all"), "Transform this selection")
      ("maxiter,m", po::value<int>(&maxiter)->default_value(5000), "Maximum number of iterations for the iterative alignment algorithm")
      ("tolerance,T", po::value<greal>(&alignment_tol)->default_value(1e-6), "Tolerance for alignment convergence")
      ("center,c", po::value<bool>(&center)->default_value(true), "Auto-center the trajectory using the alignment subset")
      ("reference,r", po::value<string>(&reference_name), "Align to a reference structure (do not use iterative method)")
      ("refsel,s", po::value<string>(&reference_sel), "Selection to align against in reference (default is the same as --align)");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Trajectory filename")
      ("prefix", po::value<string>(&prefix), "Output prefix");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);
    p.add("prefix", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj") && vm.count("prefix"))) {
      cerr << "Usage- aligner [options] model-name trajectory-name output prefix\n";
      cerr << generic;
      exit(-1);
    }

  }

  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}


void centerFrame(AtomicGroup& src, AtomicGroup& trg) {
  GCoord c = src.centroid();
  trg.translate(-c);
}

void savePDB(const string& fname, const string& meta, const AtomicGroup& grp) {
  PDB pdb = PDB::fromAtomicGroup(grp);
  pdb.renumber();
  pdb.remarks().add(meta);
  ofstream ofs(fname.c_str());
  ofs << pdb;
}




int main(int argc, char *argv[]) {

  // Parse command-line options, cache invocation header for later use...
  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  // Read the inputs...
  AtomicGroup model = createSystem(model_name);
  cout << "Read in " << model.size() << " atoms from " << model_name << endl;

  pTraj traj = createTrajectory(traj_name, model);

  cout << "Reading from trajectory " << traj_name << " with " << traj->nframes() << " frames.\n";

  // Get the selections (subsets) to operate over
  AtomicGroup align_sub = selectAtoms(model, alignment_string);

  AtomicGroup applyto_sub = selectAtoms(model, transform_string);
  applyto_sub.pruneBonds();

  cout << "Subset to align with has " << align_sub.size() << " atoms.\n";

  cout << "Subset to apply alignment transformation to has " << applyto_sub.size() << " atoms.\n";

  // Now do the alignin'...
  unsigned int nframes = traj->nframes();

  
  if (reference_name.empty()) {

    // Read in the trajectory frames and extract the coordinates for the
    // aligning subset...
    vector<AtomicGroup> frames;
    while (traj->readFrame()) {
      traj->updateGroupCoords(align_sub);
      AtomicGroup subcopy = align_sub.copy();
      frames.push_back(subcopy);
    }
    
    boost::tuple<vector<XForm>,greal, int> res = iterativeAlignment(frames, alignment_tol, maxiter);
    greal final_rmsd = boost::get<1>(res);
    cout << "Final RMSD between average structures is " << final_rmsd << endl;
    cout << "Total iters = " << boost::get<2>(res) << endl;
    
    vector<XForm> xforms = boost::get<0>(res);
    
    // Zzzzap our stored groups...
    frames.clear();
    
    cout << "Aligning transformation subset...\n";
    
    // Setup for writing DCD...
    DCDWriter dcdout(prefix + ".dcd");
    dcdout.setHeader(applyto_sub.size(), nframes, 1e-3, traj->hasPeriodicBox());
    dcdout.setTitle(header);
    dcdout.writeHeader();
    
    // Now apply the alignment transformations to the requested subsets
    for (unsigned int i = 0; i<nframes; i++) {
      traj->readFrame(i);
      traj->updateGroupCoords(model);
      applyto_sub.applyTransform(xforms[i]);
      
      if (center)
        centerFrame(align_sub, applyto_sub);
      dcdout.writeFrame(applyto_sub);
      
      if (i == 0) 
        savePDB(prefix + ".pdb", header, applyto_sub);
    }
    
  } else {
    
    AtomicGroup reference = createSystem(reference_name);
    
    if (reference_sel.empty())
      reference_sel = alignment_string;
    AtomicGroup refsub = selectAtoms(reference, reference_sel);

    if (refsub.size() != align_sub.size()) {
      cerr << boost::format("ERROR- alignment subset has %d atoms but reference subset has %d.  They must match.\n") % align_sub.size() % refsub.size();
      exit(-10);
    }

    DCDWriter dcdout(prefix + ".dcd");
    dcdout.setHeader(applyto_sub.size(), nframes, 1e-3, traj->hasPeriodicBox());
    dcdout.setTitle(header);
    dcdout.writeHeader();


    bool first = true;
    while (traj->readFrame()) {
      traj->updateGroupCoords(model);
      GMatrix M = align_sub.superposition(refsub);
      XForm W(M);
      applyto_sub.applyTransform(W);

      if (center)
        centerFrame(align_sub, applyto_sub);
      dcdout.writeFrame(applyto_sub);

      if (first) {
        savePDB(prefix + ".pdb", header, applyto_sub);
        first = false;
      }
    }

  }
}



    
