/*
  MatrixIO.hpp
  Classes for reading and writing matrices in various formats...
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



#if !defined(MATRIXWRITER_HPP)
#define MATRIXWRITER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <cassert>

#include <boost/shared_array.hpp>
#include <boost/format.hpp>

#include <StreamWrapper.hpp>
#include <Matrix.hpp>

using namespace std;

typedef unsigned int uint;


//! Class for handling writing of matrices in different formats.
template<class Derived>
class MatrixWriter {
  struct null_deleter {
    void operator()(const void *) const { }
  };

public:
  MatrixWriter(const string& s) : wfs(s) { }
  MatrixWriter(const char* s) : wfs(s) { }
  MatrixWriter(fstream& fs) : wfs(fs) { }
  
  virtual ~MatrixWriter() { }

  template<class T, class P>
  void write(const loos::Matrix<T, P>&M, const string& tag, const loos::Duple& start, const loos::Duple& end, const bool transpose = false) {
    static_cast<Derived *>(this)->write_impl(M, tag, start, end, transpose);
  }
    
  template<class T, class P>
  void write(const loos::Matrix<T, P>& M, const string& tag) {
    write(M, tag, Duple(0,0), Duple(M.rows(), M.cols()), false);
  }

  template<class T, class P>
  void write(const T* data, const string& tag, const loos::Duple& size, const loos::Duple& start, const loos::Duple& end, const bool transpose = false) {
    T* p = const_cast<T *>(data);
    boost::shared_array<T> sp(p, null_deleter());
    const Matrix<T, P> M(sp, size, size.j, size.i);
    write(M, tag, start, end, transpose);
  }

  template<class T>
  void write(const T* data, const string& tag, const int m, const int n) {
    write<T,RowMajor>(data, tag, Duple(m, n), Duple(0,0), Duple(m, n), false);
  }
  

protected:
  StreamWrapper wfs;
};
  


class RawAsciiWriter : public MatrixWriter<RawAsciiWriter> {
public:
  RawAsciiWriter(const string& s) : MatrixWriter<RawAsciiWriter>(s) { }
  RawAsciiWriter(const char* s) : MatrixWriter<RawAsciiWriter>(s) { }
  RawAsciiWriter(fstream& fs) : MatrixWriter<RawAsciiWriter>(fs) { }

  template<class T, class P>
  void write_impl(const loos::Matrix<T, P>& M, const string& tag, const loos::Duple& start, const loos::Duple& end, const bool trans) {
    *(wfs()) << boost::format("# %s\n") % tag;
    *(wfs()) << boost::format("# %d x %d (%d)\n") % M.rows() % M.cols() % trans;

    for (int j=0; j<M.rows(); j++) {
      for (int i=0; i<M.cols(); i++)
	*(wfs()) << M(j, i) << " ";
      *(wfs()) << endl;
    }
  }
};




#endif
