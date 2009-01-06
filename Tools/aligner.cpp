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

struct Globals {
  string model_name, traj_name, prefix;
  string alignment_string;
  string transform_string;

  greal alignment_tol;

  int maxiter;
  bool show_rmsd;
  bool rmsdf;
  bool center;
};


Globals globals;



void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("align,a", po::value<string>(&globals.alignment_string)->default_value("name == 'CA'"),"Align using this selection")
      ("transform,t", po::value<string>(&globals.transform_string)->default_value("all"), "Transform this selection")
      ("maxiter,m", po::value<int>(&globals.maxiter)->default_value(5000), "Maximum number of iterations for the iterative alignment algorithm")
      ("tolerance,T", po::value<greal>(&globals.alignment_tol)->default_value(1e-6), "Tolerance for alignment convergence")
      ("rmsd,r", po::value<bool>(&globals.rmsdf)->default_value(false), "Calculate RMSD against average")
      ("showrmsd,s", po::value<bool>(&globals.show_rmsd)->default_value(false), "Show RMSD for each frame")
      ("center,c", po::value<bool>(&globals.center)->default_value(true), "Auto-center the trajectory");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&globals.model_name), "Model filename")
      ("traj", po::value<string>(&globals.traj_name), "Trajectory filename")
      ("prefix", po::value<string>(&globals.prefix), "Output prefix");

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

    if (!globals.rmsdf && globals.show_rmsd) {
      globals.rmsdf = true;
      cerr << "Warning - you must use --rmsd with --showrmsd, so adding --rmsd implicitly.\n";
    }

  }

  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}



double calcRMSD(const string& octave_tag, vector<AtomicGroup>& grps) {

  unsigned int nframes = grps.size();
  double avg_rmsd = 0.0;
  if (globals.show_rmsd)
    cout << "<OCTAVE>\n";
  AtomicGroup avg = averageStructure(grps);
  if (globals.show_rmsd)
    cout << octave_tag << " = [\n";
  for (unsigned int i = 0; i<nframes; i++) {
    greal irmsd = avg.rmsd(grps[i]);
    if (globals.show_rmsd)
      cout << irmsd << " ;\n";
    avg_rmsd += irmsd;
  }
  
  if (globals.show_rmsd) {
    cout << "];\n";
    cout << "</OCTAVE>\n";
  }

  return(avg_rmsd / nframes);
}



// Group coord manipulation for calculating average structure...

void zeroCoords(AtomicGroup& g) {
  unsigned int i, n = g.size();
  
  for (i=0; i<n; i++)
    g[i]->coords() = GCoord(0,0,0);
}


void addCoords(AtomicGroup& g, const AtomicGroup& h) {
  unsigned int i, n = g.size();
  
  for (i=0; i<n; i++)
    g[i]->coords() += h[i]->coords();
}


void divCoords(AtomicGroup& g, const double d) {
  unsigned int i, n = g.size();
  
  for (i=0; i<n; i++)
    g[i]->coords() /= d;
}




void savePDB(const string& fname, const string& meta, AtomicGroup& structure) {
  PDB pdb = PDB::fromAtomicGroup(structure);
  pdb.remarks().add(meta);

  ofstream ofs(fname.c_str());
  ofs << pdb;
}




int main(int argc, char *argv[]) {

  // Parse command-line options, cache invocation header for later use...
  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  // Read the inputs...
  AtomicGroup pdb = createSystem(globals.model_name);
  cout << "Read in " << pdb.size() << " atoms from " << globals.model_name << endl;

  pTraj traj = createTrajectory(globals.traj_name, pdb);

  cout << "Reading from trajectory " << globals.traj_name << " with " << traj->nframes() << " frames.\n";

  // Get the selections (subsets) to operate over
  AtomicGroup align_sub = selectAtoms(pdb, globals.alignment_string);
  AtomicGroup applyto_sub = selectAtoms(pdb, globals.transform_string);

  align_sub.clearBonds();
  cout << "Subset to align with has " << align_sub.size() << " atoms.\n";

  applyto_sub.clearBonds();
  cout << "Subset to apply alignment transformation to has " << applyto_sub.size() << " atoms.\n";

  // Now do the alignin'...
  unsigned int nframes = traj->nframes();

  // Read in the trajectory frames and extract the coordinates for the
  // aligning subset...
  vector<AtomicGroup> frames;
  while (traj->readFrame()) {
    traj->updateGroupCoords(align_sub);
    AtomicGroup subcopy = align_sub.copy();

    if (globals.center)
      subcopy.centerAtOrigin();

    frames.push_back(subcopy);
  }

  boost::tuple<vector<XForm>,greal, int> res = iterativeAlignment(frames, globals.alignment_tol, globals.maxiter);
  greal final_rmsd = boost::get<1>(res);
  cout << "Final RMSD between average structures is " << final_rmsd << endl;
  cout << "Total iters = " << boost::get<2>(res) << endl;

  vector<XForm> xforms = boost::get<0>(res);

  if (globals.rmsdf) {
    double avg_rmsd = calcRMSD("r", frames);
    cout << "Average RMSD vs average for aligned subset = " << avg_rmsd << endl;
  }

  // Zzzzap our stored groups...
  frames.clear();

  cout << "Aligning transformation subset...\n";
  // Go ahead and make first transformation (to make VMD happy so that
  // the PDB is just a copy of the first trajectory frame...
  traj->readFrame(0);
  traj->updateGroupCoords(applyto_sub);
  AtomicGroup frame = applyto_sub.copy();
  if (globals.center)
    frame.centerAtOrigin();
  frame.applyTransform(xforms[0]);
  frame.renumber();
  savePDB(globals.prefix + ".pdb", header, frame);

  // Make the first frame...
  AtomicGroup avg = applyto_sub.copy();

  // Setup for writing DCD...
  DCDWriter dcdout(globals.prefix + ".dcd");
  dcdout.setHeader(frame.size(), nframes, 1e-3, frame.isPeriodic());
  dcdout.setTitle(header);
  dcdout.writeHeader();
  dcdout.writeFrame(frame);

  // Now apply the alignment transformations to the requested subsets
  for (unsigned int i = 1; i<nframes; i++) {
    traj->readFrame(i);
    traj->updateGroupCoords(applyto_sub);
    if (globals.center)
      applyto_sub.centerAtOrigin();

    applyto_sub.applyTransform(xforms[i]);
    frame = applyto_sub.copy();
    frame.renumber();
    dcdout.writeFrame(frame);

    addCoords(avg, applyto_sub);
  }


  divCoords(avg, nframes);

  // Write out the PDB...
  savePDB(globals.prefix + "_avg.pdb", header, avg);

  if (globals.rmsdf) {
    // Second pass to calc rmsds...

    cout << "Calculating rmsds...\n";
  
    double avg_rmsd = 0.0;
    if (globals.show_rmsd)
      cout << "<OCTAVE>\nrall = [\n";
    
    for (unsigned int i=0; i<nframes; i++) {
      traj->readFrame(i);
      traj->updateGroupCoords(applyto_sub);
      applyto_sub.applyTransform(xforms[i]);
      double rms = applyto_sub.rmsd(avg);
      if (globals.show_rmsd)
        cout << rms << " ;\n";
      avg_rmsd += rms;
    }
    
    if (globals.show_rmsd)
      cout << "];\n</OCTAVE>\n";

    avg_rmsd /= nframes;
    cout << "Average RMSD vs average for transformed subset = " << avg_rmsd << endl;
  }

}



    
