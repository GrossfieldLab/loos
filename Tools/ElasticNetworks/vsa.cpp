/*
  vsa

  (c) 2009 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

  Computes the anisotropic network model for a structure.  It does
  this by building a hessian for the structure, then computing the SVD
  of it and the corresponding pseudo-inverse (ignoring the 6 lowest
  modes).

  Usage:
    vsa environment subset radius model output_prefix

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

#include <boost/format.hpp>
#include <boost/program_options.hpp>


using namespace std;
using namespace loos;
namespace po = boost::program_options;


typedef Math::Matrix<double, Math::ColMajor> Matrix;
typedef pair<uint,uint> Range;


#if defined(__linux__)
extern "C" {
  void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
  void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
  void dggev_(char*, char*, int*, double*, int*, double*, int*, double*, double*, double*, double*, int*, double*, int*, double*, int*, int*);
}
#endif




const double normalization = 1.0;
double threshold = 1e-6;


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



Matrix submatrix(const Matrix& M, const Range& rows, const Range& cols) {
  uint m = rows.second - rows.first + 1;
  uint n = cols.second - cols.first + 1;

  Matrix A(m,n);
  for (uint i=0; i < n; ++i)
    for (uint j=0; j < m; ++j)
      A(j,i) = M(j+rows.first, i+cols.first);

  return(A);
}



Matrix mmult(Matrix& A, Matrix& B, const bool transa = false, const bool transb = false) {

  f77int m = transa ? A.cols() : A.rows();
  f77int n = transb ? B.rows() : B.cols();
  f77int k = transa ? A.rows() : A.cols();
  double alpha = 1.0;
  double beta = 0.0;

  f77int lda = transa ? k : m;
  f77int ldb = transb ? n : k;
  f77int ldc = m;


  Matrix C(m, n);

#if defined(LINUX)
  char ta = (transa ? 'T' : 'N');
  char tb = (transb ? 'T' : 'N');

  dgemm_(&ta, &tb, &m, &n, &k, &alpha, A.get(), &lda, B.get(), &ldb, &beta, C.get(), &ldc);
#else
  cblas_dgemm(CblasColMajor, transa ? CblasTrans : CblasNoTrans, transb ? CblasTrans : CblasNoTrans,
              m, n, k, alpha, A.get(), lda, B.get(), ldb, beta, C.get(), ldc);
#endif

  return(C);
}


Matrix invert(Matrix& A, const double eps = 1e-6) {

  Matrix B(A);
  boost::tuple<Matrix, Matrix, Matrix> res = svd(B);
  Matrix U(boost::get<0>(res));
  Matrix S(boost::get<1>(res));
  Matrix Vt(boost::get<2>(res));


  for (uint i=0; i<Vt.rows(); ++i) {
    double sv = S[i];
    if (sv < eps)
      sv = 0.0;
    else
      sv = 1.0 / sv;

    for (uint j=0; j<Vt.cols(); ++j)
      Vt(i,j) *= sv;
  }

  Matrix Ai = mmult(Vt, U, true, true);
  return(Ai);
}



void operator+=(Matrix& A, const Matrix& B) {
  if (A.rows() != B.rows() || A.cols() != B.cols())
    throw(logic_error("Matrices are not the same size"));

  for (uint i=0; i<A.rows() * A.cols(); ++i)
    A[i] += B[i];
}

Matrix operator+(const Matrix& A, const Matrix& B) {
  Matrix C(A);
  C += B;
  return(C);
}


void operator-=(Matrix& A, const Matrix& B) {
  if (A.rows() != B.rows() || A.cols() != B.cols())
    throw(logic_error("Matrices are not the same size"));

  for (uint i=0; i<A.rows() * A.cols(); ++i)
    A[i] -= B[i];
}

Matrix operator-(const Matrix& A, const Matrix& B) {
  Matrix C(A);
  C -= B;
  return(C);
}

Matrix eye(const uint n) {
  Matrix I(n, n);
  for (uint i=0; i<n; ++i)
    I(i,i) = 1.0;

  return(I);
}



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  uint k=1;

  string envsel(argv[k++]);
  string subsel(argv[k++]);
  double rad = strtod(argv[k++], 0);
  
  AtomicGroup model = createSystem(argv[k++]);
  string prefix(argv[k++]);

  AtomicGroup subset = selectAtoms(model, subsel);
  AtomicGroup environment = selectAtoms(model, envsel);
  AtomicGroup composite = subset + environment;

  cout << "Subset size is " << subset.size() << endl;
  cout << "Environment size is " << environment.size() << endl;

  Matrix H = hessian(composite, rad);
  cout << boost::format("Hession is %dx%d\n") % H.rows() % H.cols();
  // Now, burst out the subparts...
  uint l = subset.size() * 3;

  uint n = H.cols() - 1;
  Matrix Hss = submatrix(H, Range(0,l), Range(0,l));
  Matrix Hee = submatrix(H, Range(l+1, n), Range(l+1, n));
  Matrix Hse = submatrix(H, Range(0,l), Range(l+1, n));
  Matrix Hes = submatrix(H, Range(l+1, n), Range(0, l));

  Matrix Heei = invert(Hee);
  Matrix X = mmult(Hse, Heei);
  Matrix Y = mmult(X, Hes);
  Matrix Hssp = Hss - Y;

  Matrix Ms = eye(Hss.rows());
  Matrix Me = eye(Hee.rows());

  X = mmult(Hse, Heei);
  Y = mmult(X, Me);
  X = mmult(Y, Heei);
  Y = mmult(X, Hes);
  Matrix Msp = Ms + Y;

  // Debugging...
  cout << boost::format("Hssp is %dx%d\n") % Hssp.rows() % Hssp.cols();
  cout << boost::format("Msp is %dx%d\n") % Msp.rows() % Msp.cols();

  char jobvl = 'N';
  char jobvr = 'V';
  f77int fn = Hssp.rows();
  f77int lda = fn;
  f77int ldb = fn;
  Matrix AR(fn,1);
  Matrix AI(fn,1);
  Matrix BETA(fn,1);
  double vl;
  f77int ldvl = 1;
  Matrix VR(fn,fn);
  f77int ldvr = fn;
  f77int lwork = 64 * fn;
  Matrix WORK(lwork, 1);

  f77int info;

  cout << "Calling dggev...\n";
  dggev_(&jobvl, &jobvr, &fn, Hssp.get(), &lda, Msp.get(), &ldb, AR.get(), AI.get(), BETA.get(),
         &vl, &ldvl, VR.get(), &ldvr, WORK.get(), &lwork, &info);
  cout << "Result = " << info << endl;

  Matrix eigvals(fn, 1);
  for (int i=0; i<fn; ++i) {
    if (AI[i] == 0.0)
      eigvals[i] = AR[i] / BETA[i];
    else
      eigvals[i] = 0.0;
  }

  // normalize
  for (int i=0; i<fn; ++i) {
    double norm = 0.0;
    for (int j=0; j<fn; ++j)
      norm += (VR(j,i) * VR(j,i));
    norm = sqrt(norm);
    for (int j=0; j<fn; ++j)
      VR(j,i) /= norm;
  }

  writeAsciiMatrix(prefix + "_D.asc", eigvals, hdr);
  writeAsciiMatrix(prefix + "_U.asc", VR, hdr);
}
