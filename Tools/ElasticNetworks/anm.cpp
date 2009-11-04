/*
  anm

  (c) 2008 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

  Computes the anisotropic network model for a structure.  It does
  this by building a hessian for the structure, then computing the SVD
  of it and the corresponding pseudo-inverse (ignoring the 6 lowest
  modes).

  Usage:
    anm [selection string] radius model-name output-prefix

  Examples:
    anm 'resid >= 10 && resid <= 50 && name == "CA"' 15.0 foo.pdb foo

    This creates the following files:
      foo_H.asc      == The hessian
      foo_U.asc      == Left singular vectors
      foo_s.asc      == Singular values
      foo_V.asc      == Right singular vectors
      foo_Hi.asc     == Pseudo-inverse of H

  Notes:
    o The default selection (if none is specified) is to pick CA's
    o The output is in ASCII format suitable for use with Matlab/Octave/Gnuplot
      
  
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008,2009 Tod D. Romo
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

#include <boost/format.hpp>
#include <boost/program_options.hpp>


using namespace std;
using namespace loos;
namespace po = boost::program_options;


typedef Math::Matrix<double, Math::ColMajor> Matrix;


#if defined(__linux__)
extern "C" {
  void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
}
#endif


// Globals...  Icky poo!

// This is the Kirchoff normalization constant (see Bahar, Atilgan,
// and Erman.  Folding & Design 2:173)
double normalization = 1.0;

string selection;
string model_name;
string prefix;
double cutoff;



void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("selection,s", po::value<string>(&selection)->default_value("name == 'CA'"), "Which atoms to use for the network")
      ("cutoff,c", po::value<double>(&cutoff)->default_value(15.0), "Cutoff distance for node contact");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("prefix", po::value<string>(&prefix), "Output prefix");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("prefix", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("prefix"))) {
      cerr << "Usage- anm [options] model-name output-prefix\n";
      cerr << generic;
      exit(-1);
    }
  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}


Matrix hblock(const int i, const int j, const AtomicGroup& model, const double radius2) {

  Matrix B(3,3);
  GCoord u = model[i]->coords();
  GCoord v = model[j]->coords();
  GCoord d = v - u;

  double s = d.length2();
  if (s <= radius2) {

    for (int j=0; j<3; ++j)
      for (int i=0; i<3; ++i)
        B(i,j) = normalization * d[i]*d[j] / s;
  }

  return(B);
}



Matrix hessian(const AtomicGroup& model, const double radius) {
  
  int n = model.size();
  Matrix H(3*n,3*n);
  double r2 = radius * radius;

  for (int i=1; i<n; ++i) {
    for (int j=0; j<i; ++j) {
      Matrix B = hblock(i, j, model, r2);
      for (int x = 0; x<3; ++x)
        for (int y = 0; y<3; ++y) {
          H(i*3 + y, j*3 + x) = -B(y, x);
          H(j*3 + x, i*3 + y) = -B(x ,y);
        }
    }
  }

  // Now handle the diagonal...
  for (int i=0; i<n; ++i) {
    Matrix B(3,3);
    for (int j=0; j<n; ++j) {
      if (j == i)
        continue;
      
      for (int x=0; x<3; ++x)
        for (int y=0; y<3; ++y)
          B(y,x) += H(j*3 + y, i*3 + x);
    }

    for (int x=0; x<3; ++x)
      for (int y=0; y<3; ++y)
        H(i*3 + y, i*3 + x) = -B(y,x);
  }

  return(H);
}


int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  AtomicGroup model = createSystem(model_name);
  AtomicGroup subset = selectAtoms(model, selection);
  cerr << boost::format("Selected %d atoms from %s\n") % subset.size() % model_name;

  Timer<WallTimer> timer;
  cerr << "Calculating hessian...";
  timer.start();
  Matrix H = hessian(subset, cutoff);
  timer.stop();
  cerr << boost::format(" done (%dx%d)\n") % H.rows() % H.cols();
  cerr << timer << endl;
  writeAsciiMatrix(prefix + "_H.asc", H, header);

  // Setup for SVD calc, i.e. deal with lame f77 interface
  f77int m = subset.size() * 3;
  f77int n = m;
  f77int sn = m<n ? m : n;


  Matrix U(m,m);
  Matrix S(sn,1);
  Matrix Vt(n,n);

  char jobu = 'A';
  char jobvt = 'A';
  f77int lda = m, ldu = n, lwork = -1, info, ldvt = n;
  double prework[10], *work;

  // First, request the optimal size of the work array...
  dgesvd_(&jobu, &jobvt, &m, &n, H.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, prework, &lwork, &info);
  
  lwork = (f77int)prework[0];
  work = new double[lwork];

  // Now do the actual SVD calculation...
  cerr << "Calculating SVD - ";
  timer.start();
  dgesvd_(&jobu, &jobvt, &m, &n, H.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, work, &lwork, &info);
  timer.stop();
  cerr << "done\n";
  cerr << timer << endl;
  
  if (info > 0) {
    cerr << "Convergence error in dgesvd\n";
    exit(-3);
  } else if (info < 0) {
    cerr << "Error in " << info << "th argument to dgesvd\n";
    exit(-4);
  }

  reverseRows(S);
  reverseColumns(U);
  reverseRows(Vt);

  // Write out the parts of the SVD...
  writeAsciiMatrix(prefix + "_U.asc", U, header);
  writeAsciiMatrix(prefix + "_s.asc", S, header);
  writeAsciiMatrix(prefix + "_V.asc", Vt, header, true);

  assert(sn == n && "Insufficient singular values!");

  // Now go ahead and compute the pseudo-inverse...

  // Vt = Vt * diag(1./diag(S))
  // Remember, Vt is stored col-major but transposed, hence the
  // inverted indices...
  //
  // Note:  We have to toss the last 6 terms
  for (int i=6; i<n; i++) {
    double s = 1.0/S[i];
    for (int j=0; j<n; j++)
      Vt(i,j) *= s;
  }
  
  // Ki = Vt * U';
  // Again, Vt is internally transposed, so we have to specify
  // transposing it to dgemm in order to multiply the non-transposed
  // V...
  Matrix Hi(n,n);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, n, n, n-6, 1.0, Vt.get(), n, U.get(), n, 0.0, Hi.get(), n);
  writeAsciiMatrix(prefix + "_Hi.asc", Hi, header);
}
