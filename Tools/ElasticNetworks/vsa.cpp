/*
  vsa

  (c) 2009 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry


  Usage:
    vsa [options] subset environment model output_prefix

  See:
    Woodcock et al, J Chem Phys (2008) 129:214109
    Haffner & Zheng, J Chem Phys (2009) 130:194111

*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009 Tod D. Romo
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

#include <limits>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include "hessian.hpp"

using namespace std;
using namespace loos;
namespace po = boost::program_options;


typedef pair<uint,uint> Range;


#if defined(__linux__)
extern "C" {
  void dsygvx_(int*, char*, char*, char*, int*, double*, int*, double*, int*, double*, double*, int*, int*, double*, int*, double*, double*, int*, double*, int*, int*, int*, int*);
  void dpotrf_(char*, int*, double*, int*, int*);
  void dtrmm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);
}
#endif


// Globals...
double normalization = 1.0;
double threshold = 1e-10;

string hdr;
string subset_selection, environment_selection, model_name, prefix;
double cutoff;
int verbosity = 0;
bool debug = false;
bool occupancies_are_masses;
string psf_file;

bool parameter_free;
double power;
bool hca_method;
bool nomass;

double hca_constants[5];
bool user_defined_hca_constants(false);


void fullHelp() {

  cout <<
    "\n"
    "Computes the VSA network model given a subsystem and an\n"
    "environment selection.  The output consists of several different\n"
    "ASCII formatted matrices (that can be read by Matlab/Octave) and\n"
    "depends on whether or not masses are included in the\n"
    "calculation.  If debugging is turned on (--debug), then the\n"
    "intermediate matrices are written out:\n"
    "\tfoo_H.asc    = Composite Hessian\n"
    "\tfoo_Hss.asc  = Subsystem Hessian\n"
    "\tfoo_Hee.asc  = Environment Hessian\n"
    "\tfoo_Hse.asc  = Subsystem-Environment Hessian\n"
    "\tfoo_Heei.asc = Inverted Environment Hessian\n"
    "\tfoo_Hssp.asc = Effective Subsystem Hessian\n"
    "\tfoo_Ms.asc   = Subsystem mass (optional)\n"
    "\tfoo_Me.asc   = Environment mass (optional)\n"
    "\tfoo_Msp.asc  = Effective subsystem mass (optional)\n"
    "\tfoo_R.asc    = Cholesky decomposition of Msp (optional)\n"
    "\n\n"
    "* Unit Subsystem Mass, Zero Environment Mass *\n\n"
    "Here, the effective subsystem Hessian is created and a Singular\n"
    "Value Decomposition used to solve the eigenproblem:\n"
    "\tfoo_U.asc = Subsystem eigenvectors\n"
    "\tfoo_s.asc = Subsystem eigenvalues\n"
    "\n\n"
    "* Subsystem and Environment Mass *\n\n"
    "The generalized eigenvalue problem is solved creating the\n"
    "following matrices:\n"
    "\tfoo_Ds.asc = Subsystem eigenvalues (mass-weighted)\n"
    "\tfoo_Us.asc = Subsystem eigenvectors (mass-weighted)\n"
    "\n\n"
    "* Spring Constant Control *\n\n"
    "Different methods for assigning the spring constants in the\n"
    "Hessian can be used.  For example, \"--free 1\" selects the\n"
    "\"parameter free\" method which can be combined with the \"--power\"\n"
    "option, which controls the exponent used (the default is -2).\n"
    "Note that setting the parameter-free method does not alter the\n"
    "cutoff radius used in building the Hessian, so you may want to\n"
    "set that to something very large (i.e. \"--cutoff 100\").\n"
    "Alternatively, the \"HCA\" method can be used via the \"--hca 1\"\n"
    "option.  The constants used in HCA can be set with the\n"
    "\"--hparams r_c,k1,k2,k3,k4\" option where the spring constant, k,\n"
    "is defined as,\n"
    "\tk = k1 * s - k2        if (s <= r_c)\n"
    "\tk = k3 * pow(s, -k4)   if (s > r_c)\n"
    "and s is the distance between the nodes.\n"
    "\n\n"
    "* Mass Control *\n\n"
    "VSA, by default, assumes that masses will be present.  These can\n"
    "come from one of two sources.  If \"--psf foo.psf\" is given,\n"
    "then the masses will be assigned using the \"foo.psf\" file.  This\n"
    "assumes that the atoms are in the same order between the PSF file\n"
    "and the structure file given on the command line.  Alternatively,\n"
    "the occupancy field of the PDB can be used with the\n"
    "\"--occupancies 1\" option.  See the psf-masses tool for one way to\n"
    "copy masses into a PDB's occupancies.\n"
    "\n"
    "To disable masses (i.e. use unit masses for the subsystem and\n"
    "zero masses for the environment), use the \"--nomass 1\" option.\n"
    "\n\n"
    "* Examples *\n\n"
    "\n"
    "Compute the VSA for a transmembrane region based on segid with the\n"
    "masses stored in the occupancy field of the PDB,\n"
    "\tvsa --occupancies 1 'segid == \"TRAN\" && name == \"CA\"' 'segid != \"TRAN\" && name == \"CA\"' foo.pdb foo_vsa\n"
    "\n"
    "Compute the VSA for a transmembrane region where the selection\n"
    "is stored in a file and masses taken from a PSF file,\n"
    "\tvsa --psf foo.psf \"`cat selection` && name == 'CA'\" \"not (`cat selection`) && name == 'CA'\" foo.pdb foo_vsa\n"
    "\n"
    "Compute the mass-less VSA with CAs as the subsystem and all other\n"
    "backbone atoms as the environment,\n"
    "\tvsa --nomass 1 'name == \"CA\"' 'name =~ \"^(C|O|N)$\"' foo.pdb foo_vsa\n"
    "\n"
    "The same example as above, but using the HCA spring constants,\n"
    "\tvsa --nomass 1 --hca 1 'name == \"CA\"' 'name =~ \"^(C|O|N)$\"' foo.pdb foo_vsa\n";
    
}


void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options", 120);
    generic.add_options()
      ("help", "Produce this help message")
      ("fullhelp", "More detailed help")
      ("cutoff,c", po::value<double>(&cutoff)->default_value(15.0), "Cutoff distance for node contact")
      ("psf,p", po::value<string>(&psf_file), "Take masses from the specified PSF file")
      ("free,f", po::value<bool>(&parameter_free)->default_value(false), "Use the parameter-free method rather than a cutoff")
      ("hca,h", po::value<bool>(&hca_method)->default_value(false), "Use the HCA distance scaling method")
      ("hparams,H", po::value<string>(), "Constants to use in HCA scaling (rcut, k1, k2, k3, k4)")
      ("power,P", po::value<double>(&power)->default_value(-2.0), "Scale factor to use for parameter-free method")
      ("verbosity,v", po::value<int>(&verbosity)->default_value(0), "Verbosity level")
      ("debug,d", po::value<bool>(&debug)->default_value(false), "Turn on debugging (output intermediate matrices)")
      ("occupancies,o", po::value<bool>(&occupancies_are_masses)->default_value(false), "Atom masses are stored in the PDB occupancy field")
      ("nomass,n", po::value<bool>(&nomass)->default_value(false), "Disable mass as part of the VSA solution");


    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("subset", po::value<string>(&subset_selection), "Subset selection")
      ("env", po::value<string>(&environment_selection), "Environment selection")
      ("model", po::value<string>(&model_name), "Model filename")
      ("prefix", po::value<string>(&prefix), "Output prefix");
    
    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("subset", 1);
    p.add("env", 1);
    p.add("model", 1);
    p.add("prefix", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || vm.count("fullhelp") || !(vm.count("model") && vm.count("prefix") && vm.count("subset") && vm.count("env"))) {
      cerr << "Usage- vsa [options] subset environment model-name output-prefix\n";
      cerr << generic;
      if (vm.count("fullhelp"))
        fullHelp();
      exit(-1);
    }

    if (vm.count("hparams")) {
      string s = vm["hparams"].as<string>();
      int i = sscanf(s.c_str(), "%lf,%lf,%lf,%lf,%lf", hca_constants, hca_constants+1, hca_constants+2,
                     hca_constants+3, hca_constants+4);
      if (i != 5) {
        cerr << boost::format("Error - invalid conversion of HCA constants '%s'\n") % s;
        exit(-1);
      }
      user_defined_hca_constants = true;
    }

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}





DoubleMatrix submatrix(const DoubleMatrix& M, const Range& rows, const Range& cols) {
  uint m = rows.second - rows.first;
  uint n = cols.second - cols.first;

  DoubleMatrix A(m,n);
  for (uint i=0; i < n; ++i)
    for (uint j=0; j < m; ++j)
      A(j,i) = M(j+rows.first, i+cols.first);

  return(A);
}


boost::tuple<DoubleMatrix, DoubleMatrix> eigenDecomp(DoubleMatrix& A, DoubleMatrix& B) {

  DoubleMatrix AA = A.copy();
  DoubleMatrix BB = B.copy();

  f77int itype = 1;
  char jobz = 'V';
  char uplo = 'U';
  char range = 'I';
  f77int n = AA.rows();
  f77int lda = n;
  f77int ldb = n;
  double vl = 0.0;
  double vu = 0.0;
  f77int il = 7;
  f77int iu = n;

  char dpar = 'S';
  double abstol = 2.0 * dlamch_(&dpar);
  //double abstol = -1.0;

  f77int m;
  DoubleMatrix W(n, 1);
  DoubleMatrix Z(n, n);
  f77int ldz = n;

  f77int lwork = -1;
  f77int info;
  double *work = new double[1];

  f77int *iwork = new f77int[5*n];
  f77int *ifail = new f77int[n];

  dsygvx_(&itype, &jobz, &range, &uplo, &n, AA.get(), &lda, BB.get(), &ldb, &vl, &vu, &il, &iu, &abstol, &m, W.get(), Z.get(), &ldz, work, &lwork, iwork, ifail, &info);
  if (info != 0) {
    cerr << "ERROR- dsygvx returned " << info << endl;
    exit(-1);
  }

  lwork = work[0];
  delete[] work;
  work = new double[lwork];
  dsygvx_(&itype, &jobz, &range, &uplo, &n, AA.get(), &lda, BB.get(), &ldb, &vl, &vu, &il, &iu, &abstol, &m, W.get(), Z.get(), &ldz, work, &lwork, iwork, ifail, &info);
  if (info != 0) {
    cerr << "ERROR- dsygvx returned " << info << endl;
    exit(-1);
  }

  if (m != n-6) {
    cerr << "ERROR- only got " << m << " eigenpairs instead of " << n-6 << endl;
    exit(-10);
  }

  vector<uint> indices = sortedIndex(W);
  W = permuteRows(W, indices);
  Z = permuteColumns(Z, indices);

  boost::tuple<DoubleMatrix, DoubleMatrix> result(W, Z);
  return(result);

}


// Matrix holds column vectors.  Each vector is normalized...
void normalizeColumns(DoubleMatrix& A) {
  for (uint i=0; i<A.cols(); ++i) {
    double sum = 0.0;
    for (uint j=0; j<A.rows(); ++j)
      sum += A(j, i) * A(j, i);

    if (sum <= 0) {
      for (uint j=0; j<A.rows(); ++j)
        A(j, i) = 0.0;
    } else {
      sum = sqrt(sum);
      for (uint j=0; j<A.rows(); ++j)
        A(j, i) /= sum;
    }
  }
}



// Mass-weight eigenvectors
DoubleMatrix massWeight(DoubleMatrix& U, DoubleMatrix& M) {

  // First, compute the cholesky decomp of M (i.e. it's square-root)
  DoubleMatrix R = M.copy();
  char uplo = 'U';
  f77int n = M.rows();
  f77int lda = n;
  f77int info;
  dpotrf_(&uplo, &n, R.get(), &lda, &info);
  if (info != 0) {
    cerr << "ERROR- dpotrf() returned " << info << endl;
    exit(-1);
  }

  if (debug)
    writeAsciiMatrix(prefix + "_R.asc", R, hdr, false, ScientificMatrixFormatter<double>(24, 18));

  // Now multiply M * U
  DoubleMatrix UU = U.copy();
  f77int m = U.rows();
  n = U.cols();
  double alpha = 1.0;
  f77int ldb = m;

#if defined(__linux__)
  char side = 'L';
  char notrans = 'N';
  char diag = 'N';

  dtrmm_(&side, &uplo, &notrans, &diag, &m, &n, &alpha, R.get(), &lda, UU.get(), &ldb);

#else

  cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, alpha, R.get(), lda, UU.get(), ldb);
#endif

  normalizeColumns(UU);
  return(UU);
}
  


// Map masses from one group onto another...  Minimal error checking...
void copyMasses(AtomicGroup& target, const AtomicGroup& source) {
  if (target.size() != source.size()) {
    cerr << "ERROR- groups have different sizes in copyMasses... (maybe your PSF doesn't match the model?)\n";
    exit(-1);
  }

  for (int i=0; i<target.size(); ++i) {
    if (source[i]->name() != target[i]->name()) {
      cerr << "ERROR- atom mismatch at position " << i << endl;
      exit(-1);
    }
    target[i]->mass(source[i]->mass());
  }
}


// Copy the masses from a PSF onto a group
void massFromPSF(AtomicGroup& grp, const string& name) {
  AtomicGroup psf = createSystem(name);
  copyMasses(grp, psf);
}


// The masses are stored in the occupancy field of a PDB...
void massFromOccupancy(AtomicGroup& grp) {
  for (AtomicGroup::iterator i = grp.begin(); i != grp.end(); ++i)
      (*i)->mass((*i)->occupancy());
}


// Build the 3n x 3n diagonal mass matrix for a group
DoubleMatrix getMasses(const AtomicGroup& grp) {
  uint n = grp.size();

  DoubleMatrix M(3*n,3*n);
  for (uint i=0, k=0; i<n; ++i, k += 3) {
    M(k,k) = grp[i]->mass();
    M(k+1,k+1) = grp[i]->mass();
    M(k+2,k+2) = grp[i]->mass();
  }

  return(M);
}




int main(int argc, char *argv[]) {
  hdr = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  AtomicGroup model = createSystem(model_name);


  // Ugly way of handling multiple methods for getting masses into the equation...
  if (verbosity > 0)
    cerr << "Assigning masses...\n";

  // Get masess from somewhere, or rely on the defaults provided by
  // Atom (i.e. should be 1)
  if (! psf_file.empty())
    massFromPSF(model, psf_file);
  else if (occupancies_are_masses)
    massFromOccupancy(model);
  else if (!nomass)
    cerr << "WARNING- using default masses\n";


  // Partition the model for building the composite Hessian
  AtomicGroup subset = selectAtoms(model, subset_selection);
  AtomicGroup environment = selectAtoms(model, environment_selection);
  AtomicGroup composite = subset + environment;

  if (verbosity > 1) {
    cerr << "Subset size is " << subset.size() << endl;
    cerr << "Environment size is " << environment.size() << endl;
  }

  // Insanely high precision output
  ScientificMatrixFormatter<double> sp(24,18);

  // Determine how we're going to weight the Hessian's spring constants...
  SuperBlock* blocker = 0;
  if (parameter_free)
    blocker = new DistanceWeight(composite, power);
  else if (hca_method) {
    if (user_defined_hca_constants)
      blocker = new HCA(composite, hca_constants[0], hca_constants[1], hca_constants[2], hca_constants[3], hca_constants[4]);
    else
      blocker = new HCA(composite);
  } else
    blocker = new DistanceCutoff(composite, cutoff);

  // Now build the Hessian
  DoubleMatrix H = hessian(blocker);
  delete blocker;

  // Now, burst out the subparts...
  uint l = subset.size() * 3;
  uint n = H.cols();

  DoubleMatrix Hss = submatrix(H, Range(0,l), Range(0,l));
  DoubleMatrix Hee = submatrix(H, Range(l, n), Range(l, n));
  DoubleMatrix Hse = submatrix(H, Range(0,l), Range(l, n));
  DoubleMatrix Hes = submatrix(H, Range(l, n), Range(0, l));

  Timer<WallTimer> timer;
  if (verbosity > 0) {
    cerr << "Inverting environment hessian...\n";
    timer.start();
  }

  DoubleMatrix Heei = Math::invert(Hee);
  if (verbosity > 0) {
    timer.stop();
    cerr << timer << endl;
  }
  
  // Build the effective Hessian
  DoubleMatrix Hssp = Hss - Hse * Heei * Hes;
  if (debug) {
    writeAsciiMatrix(prefix + "_H.asc", H, hdr, false, sp);
    writeAsciiMatrix(prefix + "_Hss.asc", Hss, hdr, false, sp);
    writeAsciiMatrix(prefix + "_Hee.asc", Hee, hdr, false, sp);
    writeAsciiMatrix(prefix + "_Hse.asc", Hse, hdr, false, sp);
    writeAsciiMatrix(prefix + "_Heei.asc", Heei, hdr, false, sp);
    writeAsciiMatrix(prefix + "_Hssp.asc", Hssp, hdr, false, sp);
  }

  // Shunt in the event of using unit masses...  We can use the SVD to
  // to get the eigenpairs from Hssp
  if (nomass) {
    boost::tuple<DoubleMatrix, DoubleMatrix, DoubleMatrix> svdresult = svd(Hssp);
    DoubleMatrix U(boost::get<0>(svdresult));
    DoubleMatrix S(boost::get<1>(svdresult));

    reverseColumns(U);
    reverseRows(S);

    writeAsciiMatrix(prefix + "_U.asc", U, hdr, false, sp);
    writeAsciiMatrix(prefix + "_s.asc", S, hdr, false, sp);
    exit(0);
  }


  // Build the effective mass matrix
  DoubleMatrix Ms = getMasses(subset);
  DoubleMatrix Me = getMasses(environment);
  DoubleMatrix Msp = Ms + Hse * Heei * Me * Heei * Hes;

  if (debug) {
    writeAsciiMatrix(prefix + "_Ms.asc", Ms, hdr, false, sp);
    writeAsciiMatrix(prefix + "_Me.asc", Me, hdr, false, sp);
    writeAsciiMatrix(prefix + "_Msp.asc", Msp, hdr, false, sp);
  }

  // Run the eigen-decomposition...
  if (verbosity > 0) {
    cerr << "Running eigen-decomposition of " << Hssp.rows() << " x " << Hssp.cols() << " matrix ...";
    timer.start();
  }
  boost::tuple<DoubleMatrix, DoubleMatrix> eigenpairs;
  eigenpairs = eigenDecomp(Hssp, Msp);

  DoubleMatrix Ds = boost::get<0>(eigenpairs);
  DoubleMatrix Us = boost::get<1>(eigenpairs);

  // Need to mass-weight the eigenvectors so they're orthogonal in R3...
  if (verbosity > 0)
    cerr << "mass weighting eigenvectors...";

  DoubleMatrix MUs = massWeight(Us, Msp);

  if (verbosity > 0) {
    timer.stop();
    cerr << "done\n";
    cerr << timer << endl;
  }

  writeAsciiMatrix(prefix + "_Ds.asc", Ds, hdr, false, sp);
  writeAsciiMatrix(prefix + "_Us.asc", MUs, hdr, false, sp);
}
