/*
  average.cpp

  Computes the average structure post-aligning...
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


struct Globals {
  Globals() : align_string("name == 'CA'"),
	      avg_string("(segid != 'SOLV' && segid != 'BULK') && !hydrogen"),
	      dcdmin(0), dcdmax(0),
	      alignment_tol(1e-3)
  { }


  string align_string;
  string avg_string;
  uint dcdmin, dcdmax;
  double alignment_tol;
};



Globals globals;

static struct option long_options[] = {
  {"align", required_argument, 0, 'a'},
  {"avg", required_argument, 0, 'A'},
  {"range", required_argument, 0, 'r'},
  {0,0,0,0}
};

static const char* short_options = "a:A:r:";

void show_help(void) {
  Globals defaults;

  cout << "Usage- averager [options] <pdb> <dcd>\n";
  cout << "\t--align=string       [" << defaults.align_string << "]\n";
  cout << "\t--avg=string         [" << defaults.avg_string << "]\n";
  cout << "\t--range=min:max      [";
  if (defaults.dcdmin == 0 && defaults.dcdmax == 0)
    cout << "auto]\n";
  else
    cout << defaults.dcdmin << ":" << defaults.dcdmax << "]\n";
};


void parseOptions(int argc, char *argv[]) {
  int opt, idx;

  while ((opt = getopt_long(argc, argv, short_options, long_options, &idx)) != -1)
    switch(opt) {
    case 'A': globals.avg_string = string(optarg); break;
    case 'a': globals.align_string = string(optarg); break;
    case 'r': if (sscanf(optarg, "%u:%u", &globals.dcdmin, &globals.dcdmax) != 2) {
	cerr << "Unable to parse range.\n";
	exit(-1);
      }
      break;
    case 0: break;
    default: cerr << "Unknown option '" << (char)opt << "' - ignored.\n";
    }
}



// Group coord manipulation for calculating average structure...

void zeroCoords(AtomicGroup& g) {
  uint i, n = g.size();
  
  for (i=0; i<n; i++)
    g[i]->coords() = GCoord(0,0,0);
}


void addCoords(AtomicGroup& g, const AtomicGroup& h) {
  uint i, n = g.size();
  
  for (i=0; i<n; i++)
    g[i]->coords() += h[i]->coords();
}

void subCoords(AtomicGroup& lhs, const AtomicGroup& rhs) {
  uint i, n = lhs.size();

  for (i=0; i<n; i++)
    lhs[i]->coords() -= rhs[i]->coords();
}


void divCoords(AtomicGroup& g, const double d) {
  uint i, n = g.size();
  
  for (i=0; i<n; i++)
    g[i]->coords() /= d;
}





AtomicGroup calculateAverage(const AtomicGroup& subset, const vector<XForm>& xforms, Trajectory& dcd) {
  AtomicGroup avg = subset.copy();
  AtomicGroup frame = subset.copy();
  
  zeroCoords(avg);
  for (uint i = globals.dcdmin; i<globals.dcdmax; i++) {
    dcd.readFrame(i);
    dcd.updateGroupCoords(frame);
    frame.applyTransform(xforms[i - globals.dcdmin]);
    addCoords(avg, frame);
  }

  divCoords(avg, globals.dcdmax - globals.dcdmin);
  return(avg);
}


vector<XForm> align(const AtomicGroup& subset, Trajectory& dcd) {
  vector<AtomicGroup> frames;

  for (uint i = globals.dcdmin; i<globals.dcdmax; i++) {
    AtomicGroup frame = subset.copy();
    dcd.readFrame(i);
    dcd.updateGroupCoords(frame);
    frames.push_back(frame);
  }

  boost::tuple<vector<XForm>, greal, int> res = iterativeAlignment(frames, globals.alignment_tol, 100);
  vector<XForm> xforms = boost::get<0>(res);
  greal rmsd = boost::get<1>(res);
  int iters = boost::get<2>(res);

  cerr << "Subset alignment with " << subset.size()
       << " atoms converged to " << rmsd << " rmsd after "
       << iters << " iterations.\n";

  return(xforms);
}


int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);
  
  if (argc-optind != 2) {
    show_help();
    exit(-1);
  }

  Parser alignment_parsed(globals.align_string);
  KernelSelector align_sel(alignment_parsed.kernel());

  Parser average_parsed(globals.avg_string);
  KernelSelector avg_sel(average_parsed.kernel());

  PDB pdb(argv[optind++]);

  AtomicGroup align_subset = pdb.select(align_sel);
  if (align_subset.size() == 0) {
    cerr << "Error- no atoms selected in alignment subset.\n";
    exit(-10);
  }
  cerr << "Aligning with " << align_subset.size() << " atoms.\n";

  AtomicGroup avg_subset = pdb.select(avg_sel);
  if (avg_subset.size() == 0) {
    cerr << "Error- no atoms selected in subset to average over.\n";
    exit(-10);
  }
  cerr << "Averaging over " << avg_subset.size() << " atoms.\n";

  DCD dcd(argv[optind]);

  globals.dcdmax = (globals.dcdmax == 0) ? dcd.nframes() : globals.dcdmax+1;

  cerr << "Aligning...\n";
  vector<XForm> xforms = align(align_subset, dcd);
  cerr << "Averaging...\n";

  AtomicGroup avg = calculateAverage(avg_subset, xforms, dcd);
  
  PDB avgpdb = PDB::fromAtomicGroup(avg);
  avgpdb.remarks().add(header);
  cout << avgpdb;
}

