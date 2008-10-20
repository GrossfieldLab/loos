/*
  svd.cpp

  Computes the SVD for a trajectory.  Writes out the SVD as an
  OCTAVE-formatted text file.
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

using namespace loos;

#if defined(__linux__)
extern "C" {
  void sgesvd_(char*, char*, int*, int*, float*, int*, float*, float*, int*, float*, int*, float*, int*, int*);
}
#endif

typedef unsigned int uint;   // Ah, old-style unix C!
typedef float svdreal;

#define SVDFUNC  sgesvd_

struct Globals {
  Globals()  : alignment_string("name == 'CA'"),
	       svd_string("!(segid == 'BULK' || segid == 'SOLV')"),
	       alignment_tol(1e-6),
	       include_source(0),
	       terms(0),
	       output_prefix(""),
	       avg_name(""),
	       dcdmin(0), dcdmax(0), header("<NULL HEADER>"), mapname("") { }


  string alignment_string, svd_string;
  greal alignment_tol;
  int include_source;
  int terms;
  string output_prefix;
  string avg_name;
  uint dcdmin, dcdmax;
  string header;
  string mapname;
};




Globals globals;

static struct option long_options[] = {
  {"align", required_argument, 0, 'a'},
  {"svd", required_argument, 0, 's'},
  {"tolerance", required_argument, 0, 't'},
  {"source", no_argument, 0, 'i'},
  {"terms", required_argument, 0, 'S'},
  {"format", required_argument, 0, 'f'},
  {"prefix", required_argument, 0, 'p'},
  {"range", required_argument, 0, 'r'},
  {"avg", required_argument, 0, 'A'},
  {"map", required_argument, 0, 'm'},
  {"help", no_argument, 0, 'H'},
  {0,0,0,0}
};

static const char* short_options = "a:s:tiS:f:p:r:A:";


void show_help(void) {
  Globals defaults;

  cout << "Usage- svd [opts] pdb dcd\n";
  cout << "       --align=string       [" << defaults.alignment_string << "]\n";
  cout << "       --svd=string         [" << defaults.svd_string << "]\n";
  cout << "       --avg=fname          [" << defaults.avg_name << "]\n";
  cout << "       --tolerance=float    [" << defaults.alignment_tol << "]\n";
  cout << "       --source=bool        [" << defaults.include_source << "]\n";
  cout << "       --terms=int          [" << defaults.terms << "]\n";
  cout << "       --prefix=string      [" << defaults.output_prefix << "]\n";
  cout << "       --range=min:max      [" << defaults.dcdmin << ":" << defaults.dcdmax << "]\n";
  cout << "       --map=fname          [" << defaults.mapname << "]\n";
  cout << "       --help\n";
}


void parseOptions(int argc, char *argv[]) {
  int opt, idx;

  while ((opt = getopt_long(argc, argv, short_options, long_options, &idx)) != -1)
    switch(opt) {
    case 'A': globals.avg_name = string(optarg); break;
    case 'a': globals.alignment_string = string(optarg); break;
    case 's': globals.svd_string = string(optarg); break;
    case 't': globals.alignment_tol = strtod(optarg, 0); break;
    case 'i': globals.include_source = 1; break;
    case 'S': globals.terms = atoi(optarg); break;
    case 'H': show_help(); exit(0); break;
    case 'p': globals.output_prefix = string(optarg); break;
    case 'r': if (sscanf(optarg, "%u:%u", &globals.dcdmin, &globals.dcdmax) != 2) {
	cerr << "Unable to parse range.\n";
	exit(-1);
      }
      break;
    case 'm': globals.mapname = string(optarg); break;
    case 0: break;
    default: cerr << "Unknown option '" << (char)opt << "' - ignored.\n";
    }

}



vector<XForm> align(const AtomicGroup& subset, Trajectory& dcd) {

  boost::tuple<vector<XForm>, greal, int> res = iterativeAlignment(subset, dcd, globals.alignment_tol, 100);
  vector<XForm> xforms = boost::get<0>(res);
  greal rmsd = boost::get<1>(res);
  int iters = boost::get<2>(res);

  cerr << "Subset alignment with " << subset.size()
       << " atoms converged to " << rmsd << " rmsd after "
       << iters << " iterations.\n";

  return(xforms);
}


void writeAverage(const AtomicGroup& avg) {
  PDB avgpdb = PDB::fromAtomicGroup(avg);
  avgpdb.remarks().add(globals.header);
  ofstream ofs(globals.avg_name.c_str());
  if (!ofs)
    throw(runtime_error("Cannot open " + globals.avg_name));
  ofs << avgpdb;
}


// Calculates the transformed avg structure, then extracts the
// transformed coords from the DCD with the avg subtraced out...

float* extractCoords(const AtomicGroup& subset, const vector<XForm>& xforms, Trajectory& dcd) {
  AtomicGroup avg = averageStructure(subset, xforms, dcd);

  // Hook to get the avg structure if requested...
  if (globals.avg_name != "")
    writeAverage(avg);

  uint natoms = subset.size();
  AtomicGroup frame = subset.copy();
  uint n = globals.dcdmax - globals.dcdmin;

  float *block = new float[n * natoms * 3];
  uint cox = 0;
  for (uint i=globals.dcdmin; i<globals.dcdmax; i++) {
    dcd.readFrame(i);
    dcd.updateGroupCoords(frame);
    frame.applyTransform(xforms[i - globals.dcdmin]);

    for (uint j=0; j<natoms; j++) {
      GCoord c = frame[j]->coords() - avg[j]->coords();
      block[cox++] = c.x();
      block[cox++] = c.y();
      block[cox++] = c.z();
    }
  }

  return(block);
}



void write_map(const string& fname, const AtomicGroup& grp) {
  ofstream fout(fname.c_str());

  if (!fout) {
    cerr << "Unable to open " << fname << " for output.\n";
    exit(-10);
  }

  pAtom pa;
  AtomicGroup::Iterator iter(grp);
  int i = 0;
  while (pa = iter())
    fout << i++ << "\t" << pa->id() << "\t" << pa->resid() << endl;

}



int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);
  globals.header = header;
  parseOptions(argc, argv);

  if (argc - optind != 2) {
    cerr << "Invalid arguments.\n";
    show_help();
    exit(-1);
  }
  

  // Need to address this...
  AtomicGroup pdb = loos::createSystem(argv[optind++]);
  pTraj ptraj = loos::createTrajectory(argv[optind], pdb);
  
  // Fix max-range for DCD
  if (globals.dcdmax == 0)
    globals.dcdmax = ptraj->nframes();
  if (globals.dcdmin > ptraj->nframes() || globals.dcdmax > ptraj->nframes()) {
    cerr << "Invalid Trajectory range requested.\n";
    exit(-1);
  }

  Parser alignment_parsed(globals.alignment_string);
  KernelSelector align_sel(alignment_parsed.kernel());
  AtomicGroup alignsub = pdb.select(align_sel);
  if (alignsub.size() == 0) {
    cerr << "Error- no atoms selected to align with.\n";
    exit(-1);
  }


  Parser svd_parsed(globals.svd_string);
  KernelSelector svd_sel(svd_parsed.kernel());
  AtomicGroup svdsub = pdb.select(svd_sel);
  if (svdsub.size() == 0) {
    cerr << "Error- no atoms selected to calculate the SVD of.\n";
    exit(-1);
  }

  if (globals.mapname != "")
    write_map(globals.mapname, svdsub);

  cerr << argv[0] << ": Aligning...\n";
  vector<XForm> xforms = align(alignsub, *ptraj);
  cerr << argv[0] << ": Extracting aligned coordinates...\n";
  svdreal *A = extractCoords(svdsub, xforms, *ptraj);
  f77int m = svdsub.size() * 3;
  f77int n = globals.dcdmax - globals.dcdmin;
  f77int sn = m<n ? m : n;
  Matrix<svdreal, ColMajor> Am(A, m, n);


  if (globals.include_source)
    writeAsciiMatrix(globals.output_prefix + "A.asc", Am, header);

  double estimate = m*m*sizeof(svdreal) + n*n*sizeof(svdreal) + m*n*sizeof(svdreal) + sn*sizeof(svdreal);
  cerr << argv[0] << ": Allocating space... (" << m << "," << n << ") for " << estimate/megabytes << "Mb\n";
  char jobu = 'A', jobvt = 'A';
  f77int lda = m, ldu = m, ldvt = n, lwork= -1, info;
  svdreal prework[10], *work;

  svdreal *U = new svdreal[m*m];
  Matrix<svdreal, ColMajor> Um(U, m, m);

  svdreal *S = new svdreal[sn];
  Matrix<svdreal, ColMajor> Sm(S,sn,1);

  svdreal *Vt = new svdreal[n*n];
  Matrix<svdreal, ColMajor> Vtm(Vt,n,n);

  // First, request the optimal size of the work array...
  SVDFUNC(&jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, Vt, &ldvt, prework, &lwork, &info);
  if (info != 0) {
    cerr << "Error code from size request to dgesvd was " << info << endl;
    exit(-2);
  }

  lwork = (f77int)prework[0];
  estimate += lwork * sizeof(svdreal);
  cerr << argv[0] << ": SVD requests " << lwork << " extra space for a grand total of " << estimate / megabytes << "Mb\n";
  work = new svdreal[lwork];

  cerr << argv[0] << ": Calculating SVD...\n";
  SVDFUNC(&jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, Vt, &ldvt, work, &lwork, &info);

  if (info > 0) {
    cerr << "Convergence error in dgesvd\n";
    exit(-3);
  } else if (info < 0) {
    cerr << "Error in " << info << "th argument to dgesvd\n";
    exit(-4);
  }
  cerr << argv[0] << ": Done!\n";

  loos::writeAsciiMatrix(globals.output_prefix + "U.asc", Um, header);
  loos::writeAsciiMatrix(globals.output_prefix + "s.asc", Sm, header);
  loos::writeAsciiMatrix(globals.output_prefix + "V.asc", Vtm, header, true);

}
