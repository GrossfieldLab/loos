/*
  Matrix.hpp

  A simple matrix wrapper class.  This is not meant to be a mathematical matrix class!
*/

/*
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester
*/



#if !defined(LOOS_MATRIX_IMPL_HPP)
#define LOOS_MATRIX_IMPL_HPP

#include <iostream>
#include <ostream>
#include <string>

#include <boost/format.hpp>


// The various matrix policies require these typedefs...
typedef unsigned long ulong;
typedef unsigned int uint;


#include <MatrixOrder.hpp>           // Order (layout) policies
#include <MatrixStorage.hpp>         // Storage (physical) policies

namespace loos {
  namespace Math {

    // Forward declarations...
    template <typename T, class P, template<typename> class S>
    class Matrix;

    // Forward declarations (required for friend status) of conversion
    // functions...
    template<typename T, template<typename> class S>
    Matrix<T, RowMajor, S> reinterpretOrder(const Matrix<T, ColMajor,S>& );

    template<typename T, template<typename> class S>
    Matrix<T, ColMajor, S> reinterpretOrder(const Matrix<T, RowMajor,S>& );



    //! Templatized class for intepreting memory as a matrix stored in
    //! different formats...
    /** This class is a wrapper around a block of data allowing the user
     * to access it as a 2D matrix.  The idea is that the data can
     * internally be stored in a row-major, col-major, or triangular
     * format but the user-interface is always the same.  In addition,
     * some capacity to reinterpret the data as col or row major is
     * provided...
     *
     * Newly allocate matrices have each element initialized to 0.
     *
     * Differing physical storage implementations can be selected by
     * providing an appropriate storage policy.  The default is to use a
     * linear array in memory that is controlled by a
     * boost::shared_array pointer.  Alternatively, a SparseArray is
     * available that uses a tr1::unordered_map to implement a sparse matrix.
     *
     * Note: metadata is currently unused...
     *
     * Note: Not all functions will be available depending on how the
     * Matrix is configured via the policies.  For example, you cannot
     * call get() on a Matrix formed from a SparseArray storage policy
     * because a pointer to the map makes no sense.  On the other hand,
     * the SparseArray provides STL iterator access to the underlying
     * map which the SharedArray policy does not..
     */

    template<typename T, class OrderPolicy = ColMajor, template<typename> class StoragePolicy = SharedArray>
    class Matrix : public OrderPolicy, public StoragePolicy<T> {
    public:

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
