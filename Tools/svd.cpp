/*
  svd.cpp

  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Computes the SVD for a trajectory.  Writes out the SVD as an
  OCTAVE-formatted text file.
*/


#include <loos.hpp>
#include <pdb.hpp>
#include <Parser.hpp>
#include <Selectors.hpp>
#include <dcd.hpp>
#include <ensembles.hpp>

#include <getopt.h>
#include <cstdlib>
#include <iomanip>


typedef unsigned int uint;   // Ah, old-style unix C!

const uint kilobytes = 1024;
const uint megabytes = kilobytes * kilobytes;


struct Globals {
  Globals()  : alignment_string("name == 'CA'"),
	       svd_string("!(segid == 'BULK' || segid == 'SOLV')"),
	       alignment_tol(0.5),
	       include_source(1) { }


  string alignment_string, svd_string;
  greal alignment_tol;
  int include_source;
};




Globals globals;

static struct option long_options[] = {
  {"align", required_argument, 0, 'a'},
  {"svd", required_argument, 0, 's'},
  {"tolerance", required_argument, 0, 't'},
  {"source", no_argument, 0, 'i'},
  {0,0,0,0}
};

static const char* short_options = "a:s:ti";


void show_help(void) {
  Globals defaults;

  cout << "Usage- svd [opts] pdb dcd >output\n";
  cout << "       --align=string    [" << defaults.alignment_string << "]\n";
  cout << "       --svd=string      [" << defaults.svd_string << "]\n";
  cout << "       --tolerance=float [" << defaults.alignment_tol << "]\n";
  cout << "       --source=bool     [" << defaults.include_source << "]\n";
}


void parseOptions(int argc, char *argv[]) {
  int opt, idx;

  while ((opt = getopt_long(argc, argv, short_options, long_options, &idx)) != -1)
    switch(opt) {
    case 'a': globals.alignment_string = string(optarg); break;
    case 's': globals.svd_string = string(optarg); break;
    case 't': globals.alignment_tol = strtod(optarg, 0); break;
    case 'i': globals.include_source = 1; break;
    case 0: break;
    default: cerr << "Unknown option '" << opt << "' - ignored.\n";
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





AtomicGroup calculateAverage(const AtomicGroup& subset, const vector<XForm>& xforms, DCD& dcd) {
  uint n = dcd.nsteps();
  AtomicGroup avg = subset.copy();
  AtomicGroup frame = subset.copy();
  
  zeroCoords(avg);
  for (uint i = 0; i<n; i++) {
    dcd.readFrame(i);
    dcd.updateGroupCoords(frame);
    frame.applyTransform(xforms[i]);
    addCoords(avg, frame);
  }

  divCoords(avg, n);
  return(avg);
}


vector<XForm> align(const AtomicGroup& subset, DCD& dcd) {
  uint n = dcd.nsteps();
  vector<AtomicGroup> frames;

  for (uint i = 0; i<n; i++) {
    AtomicGroup frame = subset.copy();
    dcd.readFrame(i);
    dcd.updateGroupCoords(frame);
    frames.push_back(frame);
  }

  boost::tuple<vector<XForm>, greal, int> res = iterativeAlignment(frames, globals.alignment_tol, 100);
  vector<XForm> xforms = boost::get<0>(res);
  greal rmsd = boost::get<1>(res);
  int iters = boost::get<2>(res);

  cout << "# Subset alignment with " << subset.size()
       << " atoms converged to " << rmsd << " rmsd after "
       << iters << " iterations.\n";

  return(xforms);
}


// Calculates the transformed avg structure, then extracts the
// transformed coords from the DCD with the avg subtraced out...

double* extractCoords(const AtomicGroup& subset, const vector<XForm>& xforms, DCD& dcd) {
  AtomicGroup avg = calculateAverage(subset, xforms, dcd);
  uint n = dcd.nsteps();
  uint natoms = subset.size();
  AtomicGroup frame = subset.copy();

  double *block = new double[n * natoms * 3];
  uint cox = 0;
  for (uint i=0; i<n; i++) {
    dcd.readFrame(i);
    dcd.updateGroupCoords(frame);
    frame.applyTransform(xforms[i]);

    for (uint j=0; j<natoms; j++) {
      GCoord c = frame[j]->coords() - avg[j]->coords();
      block[cox++] = c.x();
      block[cox++] = c.y();
      block[cox++] = c.z();
    }
  }

  if (cox != n*3*natoms)
    cerr << "Ruh-Roh!\n";
  
  return(block);
}



void writeMatrix(const string& tag, const double *p, uint m, uint n, bool transpose = false) {
  uint i, j;

  cout << tag << " = [\n";
  for (j=0; j<m; j++) {
    for (i=0; i<n; i++) {
      double d;

      if (transpose)
	d = p[j*n+i];
      else
	d = p[i*m+j];

      cout << d << " ";
    }
    cout << ";\n";
  }
  cout << "];\n\n";
}



int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);
  
  cout << "# " << header << endl;
  PDB pdb(argv[optind++]);
  DCD dcd(argv[optind]);

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

  cerr << argv[0] << ": Aligning...\n";
  vector<XForm> xforms = align(alignsub, dcd);
  cerr << argv[0] << ": Extracting aligned coordinates...\n";
  double *A = extractCoords(svdsub, xforms, dcd);
  f77int m = svdsub.size() * 3;
  f77int n = dcd.nsteps();


  if (globals.include_source)
    writeMatrix("A", A, m, n);

  double estimate = (double)m*m*sizeof(double) + n*n*sizeof(double) + m*n*sizeof(double) + ((m<n)?m:n)*sizeof(double);
  cerr << argv[0] << ": Allocating space... (" << m << "," << n << ") for " << estimate/megabytes << "Mb\n";
  char jobu = 'A', jobvt = 'A';
  f77int lda = m, ldu = m, ldvt = n, lwork= -1, info;
  double prework[10], *work;

  double *U = new double[m*m];
  double *S = new double[(m<n)?m:n];
  double *Vt = new double[n*n];

  // First, request the optimal size of the work array...
  dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, Vt, &ldvt, prework, &lwork, &info);
  if (info != 0) {
    cerr << "Error code from size request to dgesvd was " << info << endl;
    exit(-2);
  }

  lwork = (f77int)prework[0];
  estimate += lwork * sizeof(double);
  cerr << argv[0] << ": SVD requests " << lwork << " extra space for a grand total of " << estimate / megabytes << "Mb\n";
  work = new double[lwork];

  cerr << argv[0] << ": Calculating SVD...\n";
  dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, Vt, &ldvt, work, &lwork, &info);

  if (info > 0) {
    cerr << "Convergence error in dgesvd\n";
    exit(-3);
  } else if (info < 0) {
    cerr << "Error in " << info << "th argument to dgesvd\n";
    exit(-4);
  }
  cerr << argv[0] << ": Done!\n";

  writeMatrix("U", U, m, m);
  writeMatrix("S", S, m, 1);
  writeMatrix("Vt", Vt, n, n, true);

  delete[] work;
  delete[] A;
  delete[] U;
  delete[] S;
  delete[] Vt;
}
