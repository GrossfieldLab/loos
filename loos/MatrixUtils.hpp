/*
  MatrixUtils.hpp

  Matrix utility functions...
*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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


#if !defined(LOOS_MATRIX_UTILS_HPP)
#define LOOS_MATRIX_UTILS_HPP

#include <string>
#include <vector>

#include <loos/loos_defs.hpp>
#include <loos/Matrix.hpp>



namespace loos {

  namespace Math {

    //! Copy one matrix into another, converting order/storage along the way.
    /**
     * Scans over all elements of a matrix.
     * This may not behave as intended with sparse matrices...
     */
    template<class T1, class P1, template<typename> class S1,
             class T2, class P2, template<typename> class S2>
    void copyMatrix(Matrix<T1,P1,S1>& A, const Matrix<T2,P2,S2>& M) {
      uint m = M.rows();
      uint n = M.cols();
      Matrix<T1,P1,S1> B(m, n);

      for (uint j=0; j<m; ++j)
        for (uint i=0; i<n; ++i)
          B(j,i) = M(j, i);

      A = B;
    }

    //! Overload for copying a sparse matrix that preserves the original
    // matrix's sparseness...
    /** Note:  static null_value SHOULD match what's in MatrixStorage, but
     *        this probably isn't guaranteed
     * Note:  Overloading templated functions is dangerous.  If this
     *        needs to be expanded, better to use an Impl class...
     */

    template<class T1, class P1, class T2, class P2>
    void copyMatrix(Matrix<T1,P1,SparseArray>& A, const Matrix<T2,P2,SparseArray>& M) {
      uint m = M.rows();
      uint n = M.cols();
      static T2 null_value;
      Matrix<T1,P1,SparseArray> B(m, n);

      for (uint j=0; j<m; ++j)
        for (uint i=0; i<n; ++i) {
          T2 t = M(j, i);
          if (t != null_value)
            B(j, i) = t;
        }
      A = B;
    }

    template<typename T, template<typename> class S>
    Matrix<T, RowMajor, S> reinterpretOrder(const Matrix<T, ColMajor, S>& A) {
      Matrix<T, RowMajor, S> result(A.rows(), A.cols());
      result.set(A);
      return(result);
    };

    template<typename T, template<typename> class S>
    Matrix<T, ColMajor, S> reinterpretOrder(const Matrix<T, RowMajor, S>& A) {
      Matrix<T, ColMajor, S> result(A.rows(), A.cols());
      result.set(A);
      return(result);
    };


    // These are left as templated functions in case we need to
    // specialize them in the future for performance reasons...

    //! Extract a row from a matrix as a vector of T
    template<typename T, class P, template<typename> class S>
    std::vector<T> getRow(const Matrix<T,P,S>& M, const uint j) {
      uint n = M.cols();
      std::vector<T> row(n);
      for (uint i=0; i<n; ++i)
        row[i] = M(j,i);

      return(row);
    }


    //! Extract a column from a matrix as a vector of T
    template<typename T, class P, template<typename> class S>
    std::vector<T> getCol(const Matrix<T,P,S>& M, const uint i) {
      uint m = M.rows();
      std::vector<T> col(m);
      for (uint j=0; j<m; ++j)
        col[j] = M(j,i);

      return(col);
    }

  }

}


#endif


