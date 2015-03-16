/*
  Matrix.hpp

  A simple matrix wrapper class.  This is not meant to be a mathematical matrix class!
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


#if !defined(LOOS_MATRIX_IMPL_HPP)
#define LOOS_MATRIX_IMPL_HPP

#include <iostream>
#include <ostream>
#include <string>

#include <loos/loos_defs.hpp>
#include <loos/boost/format.hpp>
#include <loos/MatrixOrder.hpp>           // Order (layout) policies
#include <loos/MatrixStorage.hpp>         // Storage (physical) policies

namespace loos {

  namespace Math {

    //! Specify a range for columns/rows [first,second)
    typedef std::pair<uint,uint> Range;

    // Forward declarations...
    template <typename T, class P, template<typename> class S>
    class Matrix;

    // Forward declarations (required for friend status) of conversion
    // functions...
    template<typename T, template<typename> class S>
    Matrix<T, RowMajor, S> reinterpretOrder(const Matrix<T, ColMajor,S>& );

    template<typename T, template<typename> class S>
    Matrix<T, ColMajor, S> reinterpretOrder(const Matrix<T, RowMajor,S>& );



    //! Simple matrix template class using policy classes to determine behavior
    /**
     * This class is essentially a wrapper around a block of data with
     * a matrix-style interface.  It is not (currently) meant to offer
     * the same operations that a mathematical matrix would offer,
     * despite being in the loos::Math namespace.
     *
     * There are two options to a Matrix, other than the raw data type.
     * The memory layout can be configured, i.e. row-major, column-major,
     * and triangular (symmetric).  The storage method can also be configured
     * as either a SharedArray (dense) or SparseArray (sparse).  These
     * options are passed to the Matrix template at instantiation.
     * The default options are for a column-major matrix that is dense.
     *
     * For example, the simplest declaration:
\code
loos::Math::Matrix<double> M;
\endcode
     * creates a dense, column-major matrix \a M whose elements are a double.
     *
\code
loos::Math::Matrix<float, loos::Math::RowMajor> N;
\endcode
     * This creates a dense, row-major matrix \a N whose elements are floats.
     *
     * Finally,
\code
loos::Math::Matrix<int, loos::Math::Triangular, loos::Math::SparseArray> T;
\endcode
     * creates a sparse, triangular matrix whose elements are ints.
     *
     * Access to the elements of a matrix is permitted in two different ways.
     * You can access the (j,i)'th element (j-rows, i-cols) by using operator(),
     * i.e.
\code
M(j,i) = foo;
fu = M(j+1,i);
\endcode
     * Alternatively, for dense matrices, you can access the linear array of
     * data underneath the matrix interpretation using operator[], i.e.
\code
foo = M[i*rows+j];
M[i*rows+j+1] = fu;
\endcode
     * For dense matrices, you can also access the raw block of memory by
     * getting a pointer to it,
\code
double *dblptr = M.get();
\endcode
     * 
     * Finally, you can access iterators into the matrix much like you would
     * with an STL container,
\code
std::copy(M.begin(), M.end(), ostream_iterator<double>(cout, "\n"));
\endcode
     *
     * Not all Matrix interface functions are valid for all Matrix types.
     * As mentioned above, you cannot access the pointer to a sparse matrix
     * nor use the operator[].
     *
     * Finally, all matrices are initialized with 0 for each element...
     *
     * \b NOTE
     * Although the Matrix class lives in loos::Math, the I/O
     * functions are all in the loos namespace.  The reason for this
     * is so you can import the loos namespace and get the I/O
     * functions without importing the Matrix name, thereby allowing
     * you to typedef it in your code...
\code
using namespace std;
using namespace loos;

typedef Math::Matrix<double, Math::Triangular> Matrix;
...
Matrix M;
writeAsciiMatrix(filename, M, "My matrix");
\endcode
     * 
     * For more information, see loos::writeAsciiMatrix() and loos::readAsciiMatrix()
     *
     */

    template<typename T, class OrderPolicy = ColMajor, template<typename> class StoragePolicy = SharedArray>
    class Matrix : public OrderPolicy, public StoragePolicy<T> {
    public:

      typedef T      element_type;

      //! Unitialized matrix
      Matrix() : meta("") { }

      //! Wrap an existing block of data with a Matrix.
      /**
       * This may not make sense, depending on storage policy (i.e. sparse matrices)
       */
      Matrix(T* p, const uint b, const uint a) : OrderPolicy(b, a),
                                                 StoragePolicy<T>(p, OrderPolicy::size()),
                                                 meta("") { }

      //! Create a new block of data for the requested Matrix
      Matrix(const uint b, const uint a) : OrderPolicy(b, a),
                                           StoragePolicy<T>(OrderPolicy::size()),
                                           meta("") { }

      //! Deep copy...
      Matrix<T,OrderPolicy,StoragePolicy> copy(void) const {
        Matrix<T,OrderPolicy,StoragePolicy> M(OrderPolicy::m, OrderPolicy::n);
        M.copyData(*this);

        return(M);
      }

      uint rows(void) const { return(OrderPolicy::m); }
      uint cols(void) const { return(OrderPolicy::n); }

      //! Return the appropriate element (y-rows, x-cols)
      T& operator()(const uint y, const uint x) {
        ulong i = OrderPolicy::index(y,x);
        return(StoragePolicy<T>::operator[](i));
      }

      const T& operator()(const uint y, const uint x) const {
        ulong i = OrderPolicy::index(y,x);
        return(StoragePolicy<T>::operator[](i));
      }

      void metaData(const std::string& s) { meta = s; }
      std::string metaData(void) const { return(meta); }

      //! Deallocate data...
      void reset(void) { OrderPolicy::setSize(0,0); StoragePolicy<T>::reset(); }


      //! Convert a Col-major to Row-major format
      friend Matrix<T, RowMajor, StoragePolicy> reinterpretOrder<>(const Matrix<T,ColMajor,StoragePolicy>&);

      //! Convert a Row-major to Col-major format
      friend Matrix<T, ColMajor, StoragePolicy> reinterpretOrder<>(const Matrix<T,RowMajor,StoragePolicy>&);

    private:
      std::string meta;
  
    };


    template<typename T, class P, template<typename> class S>
    std::ostream& operator<<(std::ostream& os, const Matrix<T,P,S>& M) {
    
      uint m = M.rows();
      uint n = M.cols();
      os << boost::format("# %d %d (0)\n") % m % n;
      for (uint j=0; j<m; ++j) {
        for (uint i=0; i<n; ++i)
          os << M(j, i) << " ";
        os << std::endl;
      }
      return(os);
    }

  }


}

#endif
