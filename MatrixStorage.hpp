/*
  MatrixStorage.hpp

  Storage policies for loos::Matrix
*/

/*
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester
*/



#if !defined(LOOS_MATRIX_STORAGE_HPP)
#define LOOS_MATRIX_STORAGE_HPP


#include <string>
#include <stdexcept>
#include <boost/shared_array.hpp>
#include <vector>
#include <tr1/unordered_map>


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
      SharedArray() : dim_(0), dptr(0) { }

      T* get(void) const { return(dptr.get()); }


      T& operator[](const ulong i) {
        if (i >= dim_)
          throw(std::out_of_range("Matrix index out of range"));
        return(dptr[i]);
      }

      const T& operator[](const ulong i) const {
        if (i >= dim_)
          throw(std::out_of_range("Matrix index out of range"));
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


    // There are problems with tr1::unordered_map with XCode 2.5 and
    // earlier...  We sidestep the issue by disabling.  If you've
    // built your own GCC, then chances are it's a more recent version
    // and this should work...

#if (__GNUC__ == 4 && __GNUC_MINOR__ < 1)
    // Disabled warnings...  Too much output...
    // #warning Sparse Matrix code cannot be built in this configuration
#else

    //! Storage policy for a sparse matrix (see important note in the detailed documentation).
    /**
     * This policy implements a sparse matrix via a hash.
     *
     * There are apparently some issues using the tr1::unordered_map
     * Using gcc-4.0.1.  The SparseArray policy is disabled on these
     * systems.  If you want to use it, you're going to need to build
     * a more recent version of GCC.  This issue does not affect the
     * linux systems we have tested...
     */

    template<class T>
    class SparseArray {
    public:

      typedef typename std::tr1::unordered_map<ulong, T>::const_iterator const_iterator;
      typedef typename std::tr1::unordered_map<ulong, T>::iterator iterator;

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
        typename std::tr1::unordered_map<ulong, T>::const_iterator ci = dmap.find(i);
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
      std::tr1::unordered_map<ulong, T> dmap;

    };

#endif

  }

}




#endif
