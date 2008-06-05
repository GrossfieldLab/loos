/*
  aligner.cpp

  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Aligns structures in a trajectory with a reference...

*/


#include <loos.hpp>
#include <pdb.hpp>
#include <Parser.hpp>
#include <Selectors.hpp>
#include <dcd.hpp>
#include <dcdwriter.hpp>
#include <ensembles.hpp>

#include <getopt.h>
#include <cstdlib>
#include <iomanip>

// Default values...

string alignment_string("name == 'CA'");
string transform_string("all");

greal alignment_tol = 0.5;
int maxiter = 5000;
int show_rmsd = 0;
int no_rmsd = 0;


static struct option long_options[] = {
  {"align", required_argument, 0, 'a'},
  {"transform", required_argument, 0, 't'},
  {"tolerance", required_argument, 0, 'T'},
  {"max", required_argument, 0, 'm'},
  {"show", no_argument, &show_rmsd, 1},
  {"normsd", no_argument, &no_rmsd, 1},
  {"help", no_argument, 0, 'H'},
  {0,0,0,0}
};


static const char* short_options = "a:t:T:rn";


void show_help(void) {
  cout << "Usage- aligner [opts] pdb dcd output-filename-prefix\n";
  cout << "   --align=string     [" << alignment_string << "]\n";
  cout << "   --transform=string [" << transform_string << "]\n";
  cout << "   --tolerance=float  [" << alignment_tol << "]\n";
  cout << "   --max=int          [" << maxiter << "]\n";
  cout << "   --show=bool  (rmsd)[" << show_rmsd << "]\n";
  cout << "   --normsd=boold     [" << no_rmsd << "]\n";
  cout << "   --help\n";
}


void parseOptions(int argc, char *argv[]) {
  int opt, idx;
  string al, tr;
  greal at = alignment_tol;
  int mi = maxiter;
  int sr = show_rmsd;
  int nr = no_rmsd;

  while ((opt = getopt_long(argc, argv, short_options, long_options, &idx)) != -1) {
    switch(opt) {
    case 'a': al = string(optarg); break;
    case 't': tr = string(optarg); break;
    case 'T': at = strtod(optarg, 0); break;
    case 'm': mi = atoi(optarg); break;
    case 'r': sr = 1; break;
    case 'n': nr = 1; break;
    case 'H': show_help(); exit(-1);
    case 0: break;
    default: cerr << "Unknown option '" << opt << "'\n";
    }
  }


  // Major uglage here!  Yuck!
  if (al.size() != 0)
    alignment_string = al;
  if (tr.size() != 0)
    transform_string = tr;

  alignment_tol = at;
  maxiter = mi;
  if (!show_rmsd)
    show_rmsd = sr;
  if (!no_rmsd)
    no_rmsd = nr;

}


double calcRMSD(const string& octave_tag, vector<AtomicGroup>& grps) {

  unsigned int nframes = grps.size();
  double avg_rmsd = 0.0;
  if (show_rmsd)
    cout << "<OCTAVE>\n";
  AtomicGroup avg = averageStructure(grps);
  if (show_rmsd)
    cout << octave_tag << " = [\n";
  for (unsigned int i = 0; i<nframes; i++) {
    greal irmsd = avg.rmsd(grps[i]);
    if (show_rmsd)
      cout << irmsd << " ;\n";
    avg_rmsd += irmsd;
  }
  
  if (show_rmsd) {
    cout << "];\n";
    cout << "</OCTAVE>\n";
  }

  return(avg_rmsd / nframes);
}



int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);
  if (argc - optind != 3) {
    cerr << "Invalid arguments\n";
    show_help();
    exit(-1);
  }

  // Read the inputs...
  PDB pdb(argv[optind]);
  cout << "Read in " << pdb.size() << " atoms from " << argv[optind++] << endl;
  DCD dcd(argv[optind]);
  cout << "Reading from DCD " << argv[optind++] << " with " << dcd.nsteps() << " frames.\n";
  string prefix(argv[optind]);

  // Get the selections (subsets) to operate over
  Parser alignment_parsed(alignment_string);
  KernelSelector align_sel(alignment_parsed.kernel());
  
  Parser applyto_parsed(transform_string);
  KernelSelector apply_sel(applyto_parsed.kernel());

  AtomicGroup align_sub = pdb.select(align_sel);
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
  cout << "Subset to apply alignment transformation to has " << applyto_sub.size() << " atoms.\n";

  // Now do the alignin'...
  unsigned int nframes = dcd.nsteps();

  // Read in the DCD frames and extract the coordinates for the
  // aligning subset...
  vector<AtomicGroup> frames;
  while (dcd.readFrame()) {
    dcd.updateGroupCoords(align_sub);
    AtomicGroup subcopy = align_sub.copy();
    frames.push_back(subcopy);
  }

  boost::tuple<vector<XForm>,greal, int> res = iterativeAlignment(frames, alignment_tol, maxiter);
  greal final_rmsd = boost::get<1>(res);
  cout << "Final RMSD between average structres is " << final_rmsd << endl;
  cout << "Total iters = " << boost::get<2>(res) << endl;

  vector<XForm> xforms = boost::get<0>(res);

  if (!no_rmsd) {
    double avg_rmsd = calcRMSD("r", frames);
    cout << "Average RMSD vs average for aligned subset = " << avg_rmsd << endl;
  }

  // Zzzzap our stored groups...
  frames.clear();

  // Go ahead and make first transformation (to make VMD happy so that
  // the PDB is just a copy of the first DCD frame...
  dcd.readFrame(0);
  dcd.updateGroupCoords(applyto_sub);

  // Write out the PDB...
  PDB outpdb = PDB::fromAtomicGroup(applyto_sub);
  outpdb.remarks().add(header);
  string pdb_name = prefix + ".pdb";
  ofstream ofs(pdb_name.c_str());
  ofs << outpdb;
  ofs.close();

  // Make the first frame...
  frames.push_back(applyto_sub.copy());

  // Now apply the alignment transformations to the requested subsets
  for (unsigned int i = 1; i<nframes; i++) {
    dcd.readFrame(i);
    dcd.updateGroupCoords(applyto_sub);
    applyto_sub.applyTransform(xforms[i]);
    AtomicGroup apcopy = applyto_sub.copy();    // We cache frames to
						// write everything at
						// once... 
    frames.push_back(apcopy);
  }

  DCDWriter dcdout(prefix + ".dcd", frames, outpdb.remarks().allRemarks());

  if (!no_rmsd) {
    double avg_rmsd = calcRMSD("rall", frames);
    cout << "Average RMSD vs average for transformed subset = " << avg_rmsd << endl;
  }

}



    
