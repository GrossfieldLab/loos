/*
  MatrixStorage.hpp

  Storage policies for loos::Matrix
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


#if !defined(LOOS_MATRIX_STORAGE_HPP)
#define LOOS_MATRIX_STORAGE_HPP


#include <string>
#include <stdexcept>
#include <boost/shared_array.hpp>
#include <vector>

#if __GNUC__ == 4 && __GNUC_MINOR__ < 1
#include <ext/hash_map>
#else
#include <boost/unordered_map.hpp>
#endif

#include <loos_defs.hpp>

namespace loos {
  namespace Math {

    //! Storage policy for a block of memory wrapped in a boost::shared_array pointer.
    /**
     * This policy is based on a block of memory that is wrapped in a
     * boost::shared_array pointer.  This is the policy that you need
     * to use for interfacing with Atlas.
     *
     * Handles [potentially] actual allocation of data and
     * range-checking for accesses.
     */
    template<typename T>
    class SharedArray {
    public:
      typedef const T* const_iterator;
      typedef T* iterator;

      SharedArray(const ulong n) : dim_(n) { allocate(n); }
      SharedArray(T* p, const ulong n) : dim_(n), dptr(p) { }

      // In some cases, BOOST makes dptr(0) a shared_array<int> which
      // will cause subsequent type problems.  So, we force it to be a NULL
      // pointer but with type T and wrap that...
      SharedArray() : dim_(0), dptr(static_cast<T*>(0)) { }

      T* get(void) const { return(dptr.get()); }


      T& operator[](const ulong i) {
#if defined(DEBUG)
        if (i >= dim_)
          throw(std::out_of_range("Matrix index out of range"));
#endif
        return(dptr[i]);
      }

      const T& operator[](const ulong i) const {
#if defined(DEBUG)
        if (i >= dim_)
          throw(std::out_of_range("Matrix index out of range"));
#endif
        return(dptr[i]);
      }


      iterator begin(void) { return(dptr.get()); }
      iterator end(void) { return(dptr.get() + dim_); }
    
      const_iterator begin(void) const { return(dptr.get()); }
      const_iterator end(void) const { return(dptr.get() + dim_); }



    protected:

      void set(const SharedArray<T>& s) {
        dim_ = s.dim_;
        dptr = s.dptr;
      }

      void copyData(const SharedArray<T>& s) {
        allocate(s.dim_);
        for (ulong i=0; i<dim_; ++i)
          dptr[i] = s.dptr[i];
      };

      void resize(const ulong n) {
        dim_ = n;
        allocate(n);
      }

      void reset(void) {
        dim_ = 0;
        dptr.reset();
      }


    private:

      void allocate(const ulong n) {
        dptr = boost::shared_array<T>(new T[n]);
        for (ulong i=0; i<n; ++i)
          dptr[i] = 0;
      }

      ulong dim_;
      boost::shared_array<T> dptr;

    };


    //! Storage policy for a sparse matrix (see important note in the detailed documentation).
    /**
     * This policy implements a sparse matrix via a hash.
     *
     * There are apparently some issues using the tr1::unordered_map
     * using gcc-4.0.1.  If you're using gcc-4.0.1, LOOS will switch 
     * to using the GNU provided hash_map instead.  Caveat
     * programmer...
     */

    template<class T>
    class SparseArray {
    public:

#if __GNUC__ == 4 && __GNUC_MINOR__ < 1
      typedef typename __gnu_cxx::hash_map<ulong, T>::const_iterator const_iterator;
      typedef typename __gnu_cxx::hash_map<ulong, T>::iterator iterator;
#else
      typedef typename boost::unordered_map<ulong, T>::const_iterator const_iterator;
      typedef typename boost::unordered_map<ulong, T>::iterator iterator;
#endif

      SparseArray(const ulong n) : dim_(n) { }
      SparseArray() : dim_(0) { }


      T& operator[](const ulong i) {
        if (i >= dim_)
          throw(std::out_of_range("Matrix index out of range"));
        return(dmap[i]);
      }


      // Since unordered_set::operator[] will create an entry if none
      // exists, we have to use the find m.f. to check whether or not
      // this entry has been set.  If not, then we return a const-ref to
      // a default-initialize object of type T.  This is so we can read
      // through all indices of a sparse matrix without it then
      // ballooning out to the max possible storage...

      const T& operator[](const ulong i) const {
        static T null_value;
        if (i >= dim_)
          throw(std::out_of_range("Matrix index out of range"));

#if __GNUC__ == 4 && __GNUC_MINOR__ < 1
        typename __gnu_cxx::hash_map<ulong, T>::const_iterator ci;
#else
        typename boost::unordered_map<ulong, T>::const_iterator ci;
#endif
        
        ci = dmap.find(i);
        if (ci == dmap.end()) {
          return(null_value);
        }

        return((*ci).second);
      }

      //! The actual size (# of elements) set
      ulong actualSize(void) const { return(dmap.size()); }

      // NOTE:  No get() function here since it makes no sense...

      iterator begin(void) { return(dmap.begin()); }
      iterator end(void) { return(dmap.end()); }

      const_iterator begin(void) const { return(dmap.begin()); }
      const_iterator end(void) const { return(dmap.end()); }

      //! Degree of sparseness...
      double density(void) const {
        return( (static_cast<double>(dmap.size())) / dim_ );
      }


    protected:
      void set(const SparseArray<T>& s) {
        dim_ = s.dim_;
        dmap = s.dmap;
      }

      void copyData(const SparseArray<T>& s) {
        set(s);
      }
    
      void resize(const ulong n) {
        dim_ = n;
        dmap.clear();
      };

      void reset(void) {
        resize(0);
      }

    private:
      ulong dim_;


#if __GNUC__ == 4 && __GNUC_MINOR__ < 1
      __gnu_cxx::hash_map<ulong, T> dmap;
#else
      boost::unordered_map<ulong, T> dmap;
#endif


    };
  }

}

#endif
