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
      if (info != 0)
        throw(NumericalError("Failure in SVD"));

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
      if (info != 0)
        throw(NumericalError("Failure in SVD"));


      lwork = (f77int)prework[0];
      estimate += lwork * sizeof(double);
      work = new double[lwork];

      dgesvd_(&jobu, &jobvt, &m, &n, M.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, work, &lwork, &info);

      boost::tuple<DoubleMatrix, DoubleMatrix, DoubleMatrix> result(U, S, Vt);
      return(result);
    }

    // Multiply two matrices using BLAS

    RealMatrix MMMultiply(const RealMatrix& A, const RealMatrix& B, const bool transa, const bool transb) {

      f77int m = transa ? A.cols() : A.rows();
      f77int n = transb ? B.rows() : B.cols();
      f77int k = transa ? A.rows() : A.cols();
      float alpha = 1.0;
      float beta = 0.0;

      f77int lda = transa ? k : m;
      f77int ldb = transb ? n : k;
      f77int ldc = m;


      RealMatrix C(m, n);

#if defined(__linux__)
      char ta = (transa ? 'T' : 'N');
      char tb = (transb ? 'T' : 'N');

      sgemm_(&ta, &tb, &m, &n, &k, &alpha, A.get(), &lda, B.get(), &ldb, &beta, C.get(), &ldc);
#else
      cblas_sgemm(CblasColMajor, transa ? CblasTrans : CblasNoTrans, transb ? CblasTrans : CblasNoTrans,
                  m, n, k, alpha, A.get(), lda, B.get(), ldb, beta, C.get(), ldc);
#endif

      return(C);
    }


    // Pseudo-inverse of a matrix using the SVD

    RealMatrix invert(RealMatrix& A, const float eps) {

      // The SVD (at least dgesvd) will destroy the source matrix, so we
      // need to make an explicit copy...

      RealMatrix B(A.copy());
      boost::tuple<RealMatrix, RealMatrix, RealMatrix> res = svd(B);
      RealMatrix U(boost::get<0>(res));
      RealMatrix S(boost::get<1>(res));
      RealMatrix Vt(boost::get<2>(res));


      for (uint i=0; i<Vt.rows(); ++i) {
        double sv = S[i];
        if (sv < eps)
          sv = 0.0;
        else
          sv = 1.0 / sv;

        for (uint j=0; j<Vt.cols(); ++j)
          Vt(i,j) *= sv;
      }

      RealMatrix Ai = MMMultiply(Vt, U, true, true);
      return(Ai);
    }


    void operator+=(RealMatrix& A, const RealMatrix& B) {
      if (A.rows() != B.rows() || A.cols() != B.cols())
        throw(std::logic_error("Matrices are not the same size"));

      for (uint i=0; i<A.rows() * A.cols(); ++i)
        A[i] += B[i];
    }

    RealMatrix operator+(const RealMatrix& A, const RealMatrix& B) {
      RealMatrix C(A.copy());
      C += B;
      return(C);
    }


    void operator-=(RealMatrix& A, const RealMatrix& B) {
      if (A.rows() != B.rows() || A.cols() != B.cols())
        throw(std::logic_error("Matrices are not the same size"));

      for (uint i=0; i<A.rows() * A.cols(); ++i)
        A[i] -= B[i];
    }

    RealMatrix operator-(const RealMatrix& A, const RealMatrix& B) {
      RealMatrix C(A.copy());
      C -= B;
      return(C);
    }


    void operator*=(RealMatrix& A, const float d) {
      for (uint i=0; i<A.rows() * A.cols(); ++i)
        A[i] *= d;
    }

    RealMatrix operator*(const RealMatrix& A, const float d) {
      RealMatrix B(A.copy());
      B *= d;
      return(B);
    }

    void operator*=(RealMatrix& A, const RealMatrix& B) {
      RealMatrix C = MMMultiply(A, B);
      A=C;
    }

    RealMatrix operator*(const RealMatrix& A, const RealMatrix& B) {
      RealMatrix C = MMMultiply(A, B);
      return(C);
    }


    RealMatrix operator-(RealMatrix& A) {
      RealMatrix B(A.copy());
      for (uint i=0; i<B.rows() * B.cols(); ++i)
        B[i] = -B[i];
      return(B);
    }
  

    // Create a square identity matrix

    RealMatrix eye(const uint n) {
      RealMatrix I(n, n);
      for (uint i=0; i<n; ++i)
        I(i,i) = 1.0;

      return(I);
    }


    
    RealMatrix permuteColumns(const RealMatrix& A, const std::vector<uint> indices) {
      if (indices.size() != A.cols())
        throw(std::logic_error("indices to permuteColumns must match the size of the matrix"));

      RealMatrix B(A.rows(), A.cols());

      for (uint i=0; i<A.cols(); ++i) {
        if (indices[i] > A.cols())
          throw(std::out_of_range("Permutation index is out of bounds"));
        for (uint j=0; j<A.rows(); ++j)
          B(j, i) = A(j, indices[i]);
      }

      return(B);
    }

    
    RealMatrix permuteRows(const RealMatrix& A, const std::vector<uint> indices) {
      if (indices.size() != A.rows())
        throw(std::logic_error("indices to permuteRows must match the size of the matrix"));

      RealMatrix B(A.rows(), A.cols());

      for (uint j=0; j<A.rows(); ++j) {
        if (indices[j] > A.rows())
          throw(std::out_of_range("Permutation index is out of bounds"));
        for (uint i=0; i<A.cols(); ++i)
          B(j, i) = A(indices[j], i);
      }

      return(B);
    }


    void reverseColumns(RealMatrix& A) {
      uint m = A.rows();
      uint n = A.cols();
      uint k = n / 2;

      for (uint i=0; i<k; ++i)
        for (uint j=0; j<m; ++j)
          std::swap(A(j,i), A(j,n-i-1));
    }

    void reverseRows(RealMatrix& A) {
      uint m = A.rows();
      uint n = A.cols();
      uint k = m / 2;

      for (uint j=0; j<k; ++j)
        for (uint i=0; i<n; ++i)
          std::swap(A(j, i), A(m-j-1, i));
    }


    namespace {
      double colDotProd(const RealMatrix& A, const uint i, const RealMatrix& B, const uint j) {
        double sum = 0.0;
        for (uint k=0; k<A.rows(); ++k)
          sum += A(k,i) * B(k,j);
        return(sum);
      }
    }
    
    
    double subspaceOverlap(const RealMatrix& A, const RealMatrix& B, uint nmodes) {
      if (A.rows() != B.rows())
        throw(NumericalError("Matrices have different dimensions"));

      if (nmodes == 0)
        nmodes = A.cols();
      if (nmodes > A.cols() || nmodes > B.cols())
        throw(NumericalError("Requested number of modes exceeds matrix dimensions"));
      
      double sum = 0.0;
      for (uint i=0; i<nmodes; ++i)
        for (uint j=0; j<nmodes; ++j) {
          double d = colDotProd(A, i, B, j);
          sum += d*d;
        }

      return(sum / nmodes);
    }


    double covarianceOverlap(const RealMatrix& lamA, const RealMatrix& UA, const RealMatrix& lamB, const RealMatrix& UB, const double tol) {
      if (!(UA.rows() == UB.rows() && UA.cols() == UB.cols()))
        throw(NumericalError("Matrices have different dimensions"));

      uint nmodes = UA.cols();

      double lamsum = 0.0;
      for (uint i=0; i<nmodes; ++i)
        lamsum += lamA[i] + lamB[i];

      double dblsum = 0.0;
      for (uint i=0; i<nmodes; ++i)
        for (uint j=0; j<nmodes; ++j) {
          double d = colDotProd(UA, i, UB, j);
          dblsum += sqrt(lamA[i]*lamB[j]) * d * d;
        }
      

      double num = lamsum - 2.0 * dblsum;
      double co;
      if (fabs(num) <= tol)
        co = 1.0;
      else
        co = 1.0 - sqrt( (lamsum - 2.0 * dblsum) / lamsum );

      return(co);
    }


  }

}

