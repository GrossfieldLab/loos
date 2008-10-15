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

#include <getopt.h>
#include <cstdlib>

// Default values...


struct Globals {
  Globals() : alignment_string("name == 'CA'"),
	      transform_string("name == 'CA'"),
	      alignment_tol(1e-6),
	      maxiter(5000),
	      show_rmsd(false),
	      no_rmsd(false) { }


  string alignment_string;
  string transform_string;

  greal alignment_tol;
  int maxiter;
  bool show_rmsd;
  bool no_rmsd;
};


Globals globals;




static struct option long_options[] = {
  {"align", required_argument, 0, 'a'},
  {"transform", required_argument, 0, 't'},
  {"tolerance", required_argument, 0, 'T'},
  {"max", required_argument, 0, 'm'},
  {"show", no_argument, 0, 's'},
  {"normsd", no_argument, 0, 'n'},
  {"help", no_argument, 0, 'H'},
  {0,0,0,0}
};


static const char* short_options = "a:t:T:m:snH";


void show_help(void) {
  Globals defaults;
  cout << "Usage- aligner [opts] structure-file trajectory output-filename-prefix\n";
  cout << "   --align=string     [" << defaults.alignment_string << "]\n";
  cout << "   --transform=string [" << defaults.transform_string << "]\n";
  cout << "   --tolerance=float  [" << defaults.alignment_tol << "]\n";
  cout << "   --max=int          [" << defaults.maxiter << "]\n";
  cout << "   --show       (rmsd)[" << defaults.show_rmsd << "]\n";
  cout << "   --normsd           [" << defaults.no_rmsd << "]\n";
  cout << "   --help\n";
}


void parseOptions(int argc, char *argv[]) {
  int opt, idx;

  while ((opt = getopt_long(argc, argv, short_options, long_options, &idx)) != -1) {
    switch(opt) {
    case 'a': globals.alignment_string = string(optarg); break;
    case 't': globals.transform_string = string(optarg); break;
    case 'T': globals.alignment_tol = strtod(optarg, 0); break;
    case 'm': globals.maxiter = atoi(optarg); break;
    case 's': globals.show_rmsd = true; break;
    case 'n': globals.no_rmsd = true; break;
    case 'H': show_help(); exit(-1);
    case 0: break;
    default: cerr << "Unknown option '" << (char)opt << "' - ignored.\n";
    }
  }

  if (globals.show_rmsd && globals.no_rmsd) {
    cerr << "Error- you cannot disable RMSD and then show it's results...  That's just daft!\n";
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




int main(int argc, char *argv[]) {

  // Parse command-line options, cache invocation header for later use...
  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);
  if (argc - optind != 3) {
    cerr << "Invalid arguments\n";
    show_help();
    exit(-1);
  }

  // Read the inputs...
  AtomicGroup pdb = loos::createSystem(argv[optind]);
  cout << "Read in " << pdb.size() << " atoms from " << argv[optind++] << endl;

  pTraj traj = loos::createTrajectory(argv[optind], pdb);

  cout << "Reading from trajectory " << argv[optind++] << " with " << traj->nframes() << " frames.\n";
  string prefix(argv[optind]);

  // Get the selections (subsets) to operate over
  Parser alignment_parsed(globals.alignment_string);
  KernelSelector align_sel(alignment_parsed.kernel());
  
  Parser applyto_parsed(globals.transform_string);
  KernelSelector apply_sel(applyto_parsed.kernel());

  AtomicGroup align_sub = pdb.select(align_sel);
  align_sub.clearBonds();
  if (align_sub.size() == 0) {
    cerr << "Error- there were no atoms selected from the pdb for aligning with.\n";
    exit(-10);
  }
  cout << "Subset to align with has " << align_sub.size() << " atoms.\n";

  AtomicGroup applyto_sub = pdb.select(apply_sel);
  if (applyto_sub.size() == 0) {
    cerr << "Error- there were no atoms selected from the PDB for applying the alignment to.\n";
    exit(-10);
  }
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
    frames.push_back(subcopy);
  }

  boost::tuple<vector<XForm>,greal, int> res = iterativeAlignment(frames, globals.alignment_tol, globals.maxiter);
  greal final_rmsd = boost::get<1>(res);
  cout << "Final RMSD between average structures is " << final_rmsd << endl;
  cout << "Total iters = " << boost::get<2>(res) << endl;

  vector<XForm> xforms = boost::get<0>(res);

  if (!globals.no_rmsd) {
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
  frame.applyTransform(xforms[0]);
  frame.renumber();

  // Write out the PDB...
  PDB outpdb = PDB::fromAtomicGroup(frame);
  outpdb.remarks().add(header);
  string pdb_name = prefix + ".pdb";
  ofstream ofs(pdb_name.c_str());
  ofs << outpdb;
  ofs.close();

  // Make the first frame...
  AtomicGroup avg = applyto_sub.copy();

  // Setup for writing DCD...
  DCDWriter dcdout(prefix + ".dcd");
  dcdout.setHeader(frame.size(), nframes, 1e-3, frame.isPeriodic());
  dcdout.setTitles(outpdb.remarks().allRemarks());
  dcdout.writeHeader();
  dcdout.writeFrame(frame);

  // Now apply the alignment transformations to the requested subsets
  for (unsigned int i = 1; i<nframes; i++) {
    traj->readFrame(i);
    traj->updateGroupCoords(applyto_sub);
    applyto_sub.applyTransform(xforms[i]);
    frame = applyto_sub.copy();
    frame.renumber();
    dcdout.writeFrame(frame);

    // Track average frame...
    if (!globals.no_rmsd)
      addCoords(avg, applyto_sub);
  }


  if (!globals.no_rmsd) {
    // Second pass to calc rmsds...

    cout << "Calculating rmsds...\n";
    divCoords(avg, nframes);
    
  
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



    
