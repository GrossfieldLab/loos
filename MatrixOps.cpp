/*
  MatrixOps.cpp
*/



/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo, Alan Grossfield
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



#include <MatrixOps.hpp>


namespace loos {
  namespace Math {

    boost::tuple<RealMatrix, RealMatrix, RealMatrix> svd(RealMatrix& M) {
      f77int m = M.rows();
      f77int n = M.cols();
      f77int sn = m<n ? m : n;
      long estimate = (m*m + n*n + m*n + sn) * sizeof(float);
      // This is somewhat vestigial...  We continue to estimate the
      // total storage required for the SVD in case we want to do
      // something with it in the future...



      char jobu = 'A', jobvt = 'A';
      f77int lda = m, ldu = m, ldvt = n, lwork= -1, info;
      float prework[10], *work;
    
      RealMatrix U(m, m);
      RealMatrix S(sn, 1);
      RealMatrix Vt(n, n);

      sgesvd_(&jobu, &jobvt, &m, &n, M.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, prework, &lwork, &info);
      if (info != 0) {
        // throw here...
        exit(-2);
      }

      lwork = (f77int)prework[0];
      estimate += lwork * sizeof(double);
      work = new float[lwork];

      sgesvd_(&jobu, &jobvt, &m, &n, M.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, work, &lwork, &info);

      boost::tuple<RealMatrix, RealMatrix, RealMatrix> result(U, S, Vt);
      return(result);
    }



    boost::tuple<DoubleMatrix, DoubleMatrix, DoubleMatrix> svd(DoubleMatrix& M) {
      f77int m = M.rows();
      f77int n = M.cols();
      f77int sn = m<n ? m : n;
      long estimate = (m*m + n*n + m*n + sn) * sizeof(double);


      char jobu = 'A', jobvt = 'A';
      f77int lda = m, ldu = m, ldvt = n, lwork= -1, info;
      double prework[10], *work;
    
      DoubleMatrix U(m, m);
      DoubleMatrix S(sn, 1);
      DoubleMatrix Vt(n, n);

      dgesvd_(&jobu, &jobvt, &m, &n, M.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, prework, &lwork, &info);
      if (info != 0) {
        // throw here...
        exit(-2);
      }

      lwork = (f77int)prework[0];
      estimate += lwork * sizeof(double);
      work = new double[lwork];

      dgesvd_(&jobu, &jobvt, &m, &n, M.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, work, &lwork, &info);

      boost::tuple<DoubleMatrix, DoubleMatrix, DoubleMatrix> result(U, S, Vt);
      return(result);
    }

    // Multiply two matrices using BLAS

    DoubleMatrix MMMultiply(const DoubleMatrix& A, const DoubleMatrix& B, const bool transa, const bool transb) {

      f77int m = transa ? A.cols() : A.rows();
      f77int n = transb ? B.rows() : B.cols();
      f77int k = transa ? A.rows() : A.cols();
      double alpha = 1.0;
      double beta = 0.0;

      f77int lda = transa ? k : m;
      f77int ldb = transb ? n : k;
      f77int ldc = m;


      DoubleMatrix C(m, n);

#if defined(__linux__)
      char ta = (transa ? 'T' : 'N');
      char tb = (transb ? 'T' : 'N');

      dgemm_(&ta, &tb, &m, &n, &k, &alpha, A.get(), &lda, B.get(), &ldb, &beta, C.get(), &ldc);
#else
      cblas_dgemm(CblasColMajor, transa ? CblasTrans : CblasNoTrans, transb ? CblasTrans : CblasNoTrans,
                  m, n, k, alpha, A.get(), lda, B.get(), ldb, beta, C.get(), ldc);
#endif

      return(C);
    }


    // Pseudo-inverse of a matrix using the SVD

    DoubleMatrix invert(DoubleMatrix& A, const double eps) {

      // The SVD (at least dgesvd) will destroy the source matrix, so we
      // need to make an explicit copy...

      DoubleMatrix B(A.copy());
      boost::tuple<DoubleMatrix, DoubleMatrix, DoubleMatrix> res = svd(B);
      DoubleMatrix U(boost::get<0>(res));
      DoubleMatrix S(boost::get<1>(res));
      DoubleMatrix Vt(boost::get<2>(res));


      for (uint i=0; i<Vt.rows(); ++i) {
        double sv = S[i];
        if (sv < eps)
          sv = 0.0;
        else
          sv = 1.0 / sv;

        for (uint j=0; j<Vt.cols(); ++j)
          Vt(i,j) *= sv;
      }

      DoubleMatrix Ai = MMMultiply(Vt, U, true, true);
      return(Ai);
    }


    void operator+=(DoubleMatrix& A, const DoubleMatrix& B) {
      if (A.rows() != B.rows() || A.cols() != B.cols())
        throw(std::logic_error("Matrices are not the same size"));

      for (uint i=0; i<A.rows() * A.cols(); ++i)
        A[i] += B[i];
    }

    DoubleMatrix operator+(const DoubleMatrix& A, const DoubleMatrix& B) {
      DoubleMatrix C(A.copy());
      C += B;
      return(C);
    }


    void operator-=(DoubleMatrix& A, const DoubleMatrix& B) {
      if (A.rows() != B.rows() || A.cols() != B.cols())
        throw(std::logic_error("Matrices are not the same size"));

      for (uint i=0; i<A.rows() * A.cols(); ++i)
        A[i] -= B[i];
    }

    DoubleMatrix operator-(const DoubleMatrix& A, const DoubleMatrix& B) {
      DoubleMatrix C(A.copy());
      C -= B;
      return(C);
    }


    void operator*=(DoubleMatrix& A, const double d) {
      for (uint i=0; i<A.rows() * A.cols(); ++i)
        A[i] *= d;
    }

    DoubleMatrix operator*(const DoubleMatrix& A, const double d) {
      DoubleMatrix B(A.copy());
      B *= d;
      return(B);
    }

    void operator*=(DoubleMatrix& A, const DoubleMatrix& B) {
      DoubleMatrix C = MMMultiply(A, B);
      A=C;
    }

    DoubleMatrix operator*(const DoubleMatrix& A, const DoubleMatrix& B) {
      DoubleMatrix C = MMMultiply(A, B);
      return(C);
    }


    DoubleMatrix operator-(DoubleMatrix& A) {
      DoubleMatrix B(A.copy());
      for (uint i=0; i<B.rows() * B.cols(); ++i)
        B[i] = -B[i];
      return(B);
    }
  

    // Create a square identity matrix

    DoubleMatrix eye(const uint n) {
      DoubleMatrix I(n, n);
      for (uint i=0; i<n; ++i)
        I(i,i) = 1.0;

      return(I);
    }


  }
}

