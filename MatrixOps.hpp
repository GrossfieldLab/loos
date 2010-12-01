/*
  MatrixOps.hpp
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

#if !defined MATRIXOPS_HPP
#define MATRIXOPS_HPP

#include <loos_defs.hpp>
#include <MatrixImpl.hpp>
#include <exceptions.hpp>
#include <sorting.hpp>
#include <utils.hpp>
#include <TimeSeries.hpp>
#include <stdexcept>
#include <cmath>





namespace loos {


  typedef Math::Matrix<float, Math::ColMajor> RealMatrix;
  typedef Math::Matrix<double, Math::ColMajor> DoubleMatrix;

  /**
   * Note: the operator overloads presented for DoubleMatrix are not
   * going to be efficient.  They are only provided as a convenience
   * for working with matrices in LOOS.  If you need to perform more
   * serious linear algebra operations, you are encouraged to use
   * a 3rd party library for your tool.
   */
  namespace Math {


    //! Compute the SVD of a single precision matrix
    boost::tuple<RealMatrix, RealMatrix, RealMatrix> svd(RealMatrix& M);

    //! Compute the SVD of a double precision matrix
    /**
     * The SVD functions will overwrite the source matrix \arg M
     */
    boost::tuple<DoubleMatrix, DoubleMatrix, DoubleMatrix> svd(DoubleMatrix& M);

    //! Matrix-matrix multiply (using BLAS)
    RealMatrix MMMultiply(const RealMatrix& A, const RealMatrix& B, const bool transa = false, const bool transb = false);
    DoubleMatrix MMMultiply(const DoubleMatrix& A, const DoubleMatrix& B, const bool transa = false, const bool transb = false);

    //! Pseudo-inverse of a matrix using the SVD
    RealMatrix invert(RealMatrix& A, const float eps = 1e-5);
    DoubleMatrix invert(DoubleMatrix& A, const double eps = 1e-5);

    //! An identity matrix of size n
    RealMatrix eye(const uint n);
    DoubleMatrix deye(const uint n);


    //! Extracts a submatrix
    loos::DoubleMatrix submatrix(const loos::DoubleMatrix& M, const Range& rows, const Range& cols);

    //! Normalizes each column as a column-vector
    void normalizeColumns(loos::DoubleMatrix& A);

    //! Overloaded operators for RealMatrix matrices (see important note below)
    void operator+=(RealMatrix& A, const RealMatrix& B);
    RealMatrix operator+(const RealMatrix& A, const RealMatrix& B);
    void operator-=(RealMatrix& A, const RealMatrix& B);
    RealMatrix operator-(const RealMatrix& A, const RealMatrix& B);
    void operator*=(RealMatrix& A, const float d);
    RealMatrix operator*(const RealMatrix& A, const float d);
    void operator*=(RealMatrix& A, const RealMatrix& B);
    RealMatrix operator*(const RealMatrix& A, const RealMatrix& B);
    RealMatrix operator-(RealMatrix& A);

    void operator+=(DoubleMatrix& A, const DoubleMatrix& B);
    DoubleMatrix operator+(const DoubleMatrix& A, const DoubleMatrix& B);
    void operator-=(DoubleMatrix& A, const DoubleMatrix& B);
    DoubleMatrix operator-(const DoubleMatrix& A, const DoubleMatrix& B);
    void operator*=(DoubleMatrix& A, const double d);
    DoubleMatrix operator*(const DoubleMatrix& A, const double d);
    void operator*=(DoubleMatrix& A, const DoubleMatrix& B);
    DoubleMatrix operator*(const DoubleMatrix& A, const DoubleMatrix& B);
    DoubleMatrix operator-(DoubleMatrix& A);

    //! Returns a copy of the matrix with the columns permuted by the indices
    template<typename T>
    T permuteColumns(const T& A, const std::vector<uint>& indices) {
      if (indices.size() != A.cols())
        throw(std::logic_error("indices to permuteColumns must match the size of the matrix"));

      T B(A.rows(), A.cols());

      for (uint i=0; i<A.cols(); ++i) {
        if (indices[i] > A.cols())
          throw(std::out_of_range("Permutation index is out of bounds"));
        for (uint j=0; j<A.rows(); ++j)
          B(j, i) = A(j, indices[i]);
      }

      return(B);
    }

    //! Returns a copy of the matrix with the rows permuted by the indices
    template<typename T>
    T permuteRows(const T& A, const std::vector<uint>& indices) {
      if (indices.size() != A.rows())
        throw(std::logic_error("indices to permuteRows must match the size of the matrix"));

      T B(A.rows(), A.cols());

      for (uint j=0; j<A.rows(); ++j) {
        if (indices[j] > A.rows())
          throw(std::out_of_range("Permutation index is out of bounds"));
        for (uint i=0; i<A.cols(); ++i)
          B(j, i) = A(indices[j], i);
      }

      return(B);
    }

    
    //! Randomly shuffle the columns of a matrix
    template<typename T>
    T shuffleColumns(const T& A) {
      std::vector<float> random_numbers(A.cols());
      base_generator_type& rng = rng_singleton();
      boost::uniform_real<> rngmap(0.0, 1.0);
      boost::variate_generator<base_generator_type&, boost::uniform_real<> > rnd(rng, rngmap);

      for (uint i=0; i<A.cols(); ++i)
        random_numbers[i] = rnd();

      std::vector<uint> indices = sortedIndex(random_numbers);
      return(permuteColumns(A, indices));
    }


    //!! Randomly shuffle the rows of a single column vector
    template<typename T>
    T shuffleColumnVector(const T& v) {
      std::vector<float> random_numbers(v.size());
      base_generator_type& rng = rng_singleton();
      boost::uniform_real<> rngmap(0.0, 1.0);
      boost::variate_generator<base_generator_type&, boost::uniform_real<> > rnd(rng, rngmap);

      for (uint i=0; i<v.size(); ++i)
        random_numbers[i] = rnd();

      std::vector<uint> indices = sortedIndex(random_numbers);

      // This forces the result to -always- be a column vector...
      T u(v.size(), 1);
      for (uint i=0; i<v.size(); ++i)
        u[i] = v[indices[i]];

      return(u);
    }


    template<typename T>
    void reverseColumns(T& A) {
      uint m = A.rows();
      uint n = A.cols();
      uint k = n / 2;

      for (uint i=0; i<k; ++i)
        for (uint j=0; j<m; ++j)
          std::swap(A(j,i), A(j,n-i-1));
    }

    template<typename T>
    void reverseRows(T& A) {
      uint m = A.rows();
      uint n = A.cols();
      uint k = m / 2;

      for (uint j=0; j<k; ++j)
        for (uint i=0; i<n; ++i)
          std::swap(A(j, i), A(m-j-1, i));
    }



    namespace {
      template<typename T>
      double colDotProd(const T& A, const uint i, const T& B, const uint j) {
        double sum = 0.0;
        for (uint k=0; k<A.rows(); ++k)
          sum += A(k,i) * B(k,j);
        return(sum);
      }
    }
    
    template<typename T>
    double subspaceOverlap(const T& A, const T& B, uint nmodes = 0) {
      if (A.rows() != B.rows())
        throw(NumericalError("subspaceOverlap: Matrices have different dimensions"));

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



    //! Computes the covariance overlap between two subspaces
    /**
     * This function expects a set of eigenpairs for comparison.  The
     * eigenvalues are stored in a Matrix-vector (i.e. an nx1
     * matrix).  The eigenvectors are stored in the columns of the
     * respective matrices.
     *
     * Note: Be sure that the eigenvalues are scaled appropriately.
     * For example, when comparing PCA and ENM, the ENM eigenvalues
     * are inversely proportional to the PCA eigenvalues while the PCA
     * "eigenvalues" are actually the singular values (and hence the
     * square-root of the eigenvalues of AA'
     *
     * Note: It is possible for double sum to be slightly greater than
     * 2x the sum of the eigenvalues, which results in trying to take
     * the square root of a negative number.  To prevent this, we
     * actually use the absolute value of the difference.
     *
     * Note: Due to rounding errors in single precision, it is
     * possible that the covariance overlap of a set of eigenpairs
     * against itself will not come out to be exactly 1, but will be
     * close (i.e. to within 1e-3).
     */
    template<typename T>
    double covarianceOverlap(const T& lamA, const T& UA, const T& lamB, const T& UB) {
      if (!(UA.rows() == UB.rows() && lamA.rows() <= UA.cols() && lamB.rows() <= UB.cols()))
        throw(NumericalError("covarianceOverlap: Matrices have incorrect dimensions"));

      // X = abs(UB'*UA)
      T X = MMMultiply(UB,UA,true,false);
      for (ulong i = 0; i < X.size(); ++i)
        X[i] = fabs(X[i]);
      
      // L = abs(lamB*lamA')
      T L = MMMultiply(lamB, lamA, false, true);

      // y = sum(sum(sqrt(L).*X.*X));
      double y = 0;
      for (ulong i = 0; i<X.size(); ++i)
        y += sqrt(L[i]) * X[i] * X[i];

      // e = sum(s.*t);
      double e =0;
      for (ulong i=0; i<lamA.size(); ++i)
        e += lamA[i] + lamB[i];

      double num = e - 2.0 * y;
      double co = 1.0 - sqrt( fabs(num) / e );

      return(co);
    }


    template<typename T>
    boost::tuple<double, double> zCovarianceOverlap(const T& lamA, const T& UA, const T& lamB, const T& UB, const uint tries) {
      double coverlap = covarianceOverlap(lamA, UA, lamB, UB);
      std::vector<double> random_coverlaps(tries);

      for (uint i=0; i<tries; ++i) {
        T shuffled_lamA = shuffleColumnVector(lamA);
        T shuffled_lamB = shuffleColumnVector(lamB);
        random_coverlaps[i] = covarianceOverlap(shuffled_lamA, UA, shuffled_lamB, UB);
      }

      TimeSeries<double> ts(random_coverlaps);
      double score = (coverlap - ts.average()) / ts.stdev();

      boost::tuple<double, double> result(score, coverlap);
      return(result);
    }


  };


};


#endif
