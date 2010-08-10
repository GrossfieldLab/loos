//
// vsa-lib
// (c) 2010 Tod D. Romo
// 

#include "vsa-lib.hpp"

using namespace std;
using namespace loos;


boost::tuple<DoubleMatrix, DoubleMatrix> VSA::eigenDecomp(DoubleMatrix& A, DoubleMatrix& B) {

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



// Mass-weight eigenvectors
DoubleMatrix VSA::massWeight(DoubleMatrix& U, DoubleMatrix& M) {

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

  if (debugging_)
    writeAsciiMatrix(prefix_ + "_R.asc", R, meta_, false);

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



void VSA::solve() {
  buildHessian();
    
  uint n = hessian_.cols();
  uint l = subset_size_ * 3;
    
  DoubleMatrix Hss = submatrix(hessian_, Range(0,l), Range(0,l));
  DoubleMatrix Hee = submatrix(hessian_, Range(l, n), Range(l, n));
  DoubleMatrix Hse = submatrix(hessian_, Range(0,l), Range(l, n));
  DoubleMatrix Hes = submatrix(hessian_, Range(l, n), Range(0, l));

  if (debugging_) {
    writeAsciiMatrix(prefix_ + "_H.asc", hessian_, meta_, false);
    writeAsciiMatrix(prefix_ + "_Hss.asc", Hss, meta_, false);
    writeAsciiMatrix(prefix_ + "_Hee.asc", Hee, meta_, false);
    writeAsciiMatrix(prefix_ + "_Hse.asc", Hse, meta_, false);
  }

  DoubleMatrix Heei = Math::invert(Hee);
  
  // Build the effective Hessian
  Hssp_ = Hss - Hse * Heei * Hes;

  if (debugging_)
    writeAsciiMatrix(prefix_ + "_Hssp.asc", Hssp_, meta_, false);

  // Shunt in the event of using unit masses...  We can use the SVD to
  // to get the eigenpairs from Hssp
  if (masses_.rows() == 0) {
    boost::tuple<DoubleMatrix, DoubleMatrix, DoubleMatrix> svdresult = svd(Hssp_);

    eigenvecs_ = boost::get<0>(svdresult);
    eigenvals_ = boost::get<1>(svdresult);

    reverseColumns(eigenvecs_);
    reverseRows(eigenvals_);
    return;

  }


  // Build the effective mass matrix
  DoubleMatrix Ms = submatrix(masses_, Range(0, l), Range(0, l));
  DoubleMatrix Me = submatrix(masses_, Range(l, n), Range(l, n));

  Msp_ = Ms + Hse * Heei * Me * Heei * Hes;

  if (debugging_) {
    writeAsciiMatrix(prefix_ + "_Ms.asc", Ms, meta_, false);
    writeAsciiMatrix(prefix_ + "_Me.asc", Me, meta_, false);
    writeAsciiMatrix(prefix_ + "_Msp.asc", Msp_, meta_, false);
  }

  // Run the eigen-decomposition...
  boost::tuple<DoubleMatrix, DoubleMatrix> eigenpairs;
  eigenpairs = eigenDecomp(Hssp_, Msp_);

  eigenvals_ = boost::get<0>(eigenpairs);
  DoubleMatrix Us = boost::get<1>(eigenpairs);

  // Need to mass-weight the eigenvectors so they're orthogonal in R3...
  eigenvecs_ = massWeight(Us, Msp_);
}
