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

    //! Pseudo-inverse of a matrix using the SVD
    RealMatrix invert(RealMatrix& A, const float eps = 1e-5);

    //! An identity matrix of size n
    RealMatrix eye(const uint n);


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

    //! Returns a copy of the matrix with the columns permuted by the indices
    RealMatrix permuteColumns(const RealMatrix& A, const std::vector<uint> indices);
    //! Returns a copy of the matrix with the rows permuted by the indices
    RealMatrix permuteRows(const RealMatrix& A, const std::vector<uint> indices);
    
    //! Reverses the columns in place
    void reverseColumns(RealMatrix& A);
    //! Reverses the rows in place
    void reverseRows(RealMatrix& A);

    //! Computes the overlap between two subspaces (matrices of column vectors)
    double subspaceOverlap(const RealMatrix& A, const RealMatrix& B, uint nmodes = 0);

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
     */
    double covarianceOverlap(const RealMatrix& lamA, const RealMatrix& UA, const RealMatrix& lamB, const RealMatrix& UB, const double tol = 1e-3, uint nmodes = 0);
  };


};


#endif
