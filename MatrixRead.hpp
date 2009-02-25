/*
  MatrixRead.hpp

  Reading of Matrix objects...
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



#if !defined(MATRIXREAD_HPP)
#define MATRIXREAD_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>   // ?
#include <stdexcept>
#include <cassert>
#include <iterator>

#include <utility>

#include <boost/format.hpp>

#include <Matrix.hpp>


namespace loos {

  //! Generic reading error class
  class MatrixReadError : public std::runtime_error {
  public:
    explicit MatrixReadError(const std::string& msg) : runtime_error(msg) { }
  };

  // Forward declaration for reading implementations...
  template<class T, class P, template<typename> class S>
  struct MatrixReadImpl;



  // The following are the templated global functions.  Do not
  // overload/specialize them.  Instead, specialize the
  // MatrixReadImpl class...


  //! Read in a matrix from a stream returning a newly created matrix
  template<class T, class P, template<typename> class S>
  Math::Matrix<T,P,S> readAsciiMatrix(std::istream& is) {
    return(MatrixReadImpl<T,P,S>::read(is));
  }

  //! Read in a matrix from a stream storing it in the specified matrix
  template<class T, class P, template<typename> class S>
  void readAsciiMatrix(std::istream& is, Math::Matrix<T,P,S>& M) {
    M = MatrixReadImpl<T,P,S>::read(is);
  }

  //! Read in a matrix from a file returning a newly created matrix
  template<class T, class P, template<typename> class S>
  Math::Matrix<T,P,S> readAsciiMatrix(const std::string& fname) {
    std::ifstream ifs(fname.c_str());
    if (!ifs)
      throw(MatrixReadError("Cannot open " + fname + " for reading."));
    return(MatrixReadImpl<T,P,S>::read(ifs));
  }

  //! Read in a matrix from a file storing it in the specified matrix
  template<class T, class P, template<typename> class S>
  void readAsciiMatrix(const std::string& fname, Math::Matrix<T,P,S>& M) {
    std::ifstream ifs(fname.c_str());
    if (!ifs)
      throw(MatrixReadError("Cannot open " + fname + " for reading."));
    M = MatrixReadImpl<T,P,S>::read(ifs);
  }

  // Implementations and specializations...

  template<class T, class P, template<typename> class S>
  struct MatrixReadImpl {
    static Math::Matrix<T,P,S> read(std::istream& is) {
      std::string inbuf;
      int m = 0, n = 0;

      // First, search for the marker...
      int k = 0;
      while (getline(is, inbuf).good()) {
        ++k;
        int i = sscanf(inbuf.c_str(), "# %d %d", &m, &n);
        if (i == 2)
          break;
      }
      if (m == 0 && n == 0)
        throw(MatrixReadError("Could not find magic marker in matrix file"));
      if (m == 0 || n == 0)
        throw(MatrixReadError("Error while reading magic marker"));


      T datum;
      Math::Matrix<T,P,S> R(m, n);
      for (int j=0;j < m; j++) {
        for (int i=0; i<n; i++) {
          if (!(is >> datum)) {
            std::stringstream s;
            s << "Invalid conversion on matrix read at (" << j << "," << i << ")";
            throw(MatrixReadError(s.str()));
          }
          R(j, i) = datum;
        }
      }

      return(R);
    }
  };

  //! Special handling for sparse matrices
  template<class T, class P>
  struct MatrixReadImpl<T,P,Math::SparseArray> {
    static Math::Matrix<T,P, Math::SparseArray> read(std::istream& is) {
      std::string inbuf;
      uint m = 0, n = 0;
      ulong l = 0;

      // First, search for the marker...
      while (getline(is, inbuf).good()) {
        char buf[20];
        int i = sscanf(inbuf.c_str(), "# %lu %u %u %10s", &l, &m, &n, buf);
        if (i != 4)
          continue;
        if (strncmp(buf, "SPARSE", 6) == 0)
          break;
        throw(MatrixReadError("Magic matrix line found, but the matrix appears not to be sparse."));
      }

      if (m == 0)
        throw(MatrixReadError("Could not find magic matrix line"));
    
      Math::Matrix<T, P, Math::SparseArray> M(m, n);
      for (ulong i = 0; i<l; ++i) {
        ulong j;
        T datum;

        // Get the index
        if (!(is >> j)) {
          std::stringstream s;
          s << "Invalid conversion on matrix read at [" << i << "]";
          throw(MatrixReadError(s.str()));
        }

        // Get the value
        if (!(is >> datum)) {
          std::stringstream s;
          s << "Invalid conversion on matrix read at [" << i << "]";
          throw(MatrixReadError(s.str()));
        }
        M[j] = datum;
      }

      return(M);
    }
  };


  //! Special handling for triangular matrices
  template<class T, template<typename> class S>
  struct MatrixReadImpl<T,Math::Triangular,S> {
    static Math::Matrix<T, Math::Triangular, S> read(std::istream& is) {
      std::string inbuf;
      int m = 0;

      // First, search for the marker...
      while (getline(is, inbuf).good()) {
        char buf[20];
        int i = sscanf(inbuf.c_str(), "# %d %10s", &m, buf);
        if (i != 2)
          continue;
        if (strncmp(buf, "TRIANGULAR", 10) == 0)
          break;
        throw(MatrixReadError("Magic matrix line found, but the matrix appears not to be triangular."));
      }

      if (m == 0)
        throw(MatrixReadError("Could not find magic matrix line"));
      int n = m;

      Math::Matrix<T, Math::Triangular,S> R(m, n);
      long s = R.size();
      T datum;

      for (long i=0; i<s; i++) {
        if (!(is >> datum)) {
          std::stringstream s;
          s << "Invalid conversion on matrix read at [" << i << "]";
          throw(MatrixReadError(s.str()));
        }
        R[i] = datum;
      }

      return(R);
    }


  };
  

}

#endif
