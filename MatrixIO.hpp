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


#if !defined(MATRIXIO_HPP)
#define MATRIXIO_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <cassert>

#include <boost/shared_array.hpp>
#include <boost/format.hpp>

#include <Matrix.hpp>


namespace loos {

  // --- Writing ---

  template<class T, class P>
  ostream& writeAsciiMatrix(ostream& os, const Matrix<T,P>& M, const string& meta, const Duple& start, const Duple& end, const bool trans = false) {
    os << "# " << meta << endl;
    os << boost::format("# %d x %d (%d)\n") % M.rows() % M.cols() % trans;
    for (int j=start.j; j<end.j; j++) {
      for (int i=start.i; i<end.i; i++)
	if (trans)
	  os << M(i, j) << " ";
	else
	  os << M(j, i) << " ";
      os << endl;
    }
    return(os);
  }

  template<class T, class P>
  ostream& writeAsciiMatrix(ostream& os, const Matrix<T, P>& M, const string& meta, const bool trans = false) {
    Duple start(0,0);
    Duple end(M.rows(), M.cols());
    writeAsciiMatrix(os, M, meta, start, end, trans);
    return(os);
  }

  template<class T, class P>
  void writeAsciiMatrix(const string& fname, const Matrix<T,P>& M, const string& meta, const Duple& start, const Duple& end, const bool trans = false) {
    ofstream ofs(fname.c_str());
    if (ofs == 0)
      throw(runtime_error("Cannot open " + fname + " for writing."));
    writeAsciiMatrix(ofs, M, meta, start, end, trans);
  }

  template<class T, class P>
  void writeAsciiMatrix(const string& fname, const Matrix<T,P>& M, const string& meta, const bool trans = false) {
    Duple start(0,0);
    Duple end(M.rows(), M.cols());
    loos::writeAsciiMatrix(fname, M, meta, start, end, trans);
  }

  // --- Specializations for writing ---

  template<class T>
  ostream& writeAsciiMatrix(ostream& os, const Matrix<T, Triangular>& M, const string& meta) {
    os << "# " << meta << endl;
    os << boost::format("# %d TRIANGULAR\n") % M.rows();
    long s = M.size();
    for (long i=0; i<s; i++)
      os << M[i] << endl;

    return(os);
  }

  template<class T>
  void writeAsciiMatrix(const string& fname, const Matrix<T, Triangular>& M, const string& meta) {
    ofstream ofs(fname.c_str());
    if (ofs == 0)
      throw(runtime_error("Cannot open " + fname + " for writing."));
    writeAsciiMatrix(ofs, M, meta);
  }

  template<class T>
  ostream& writeTriangularAsFull(ostream& os, const Matrix<T, Triangular>& M, const string& meta) {
    int n = M.rows();
    os << "# " << meta << endl;
    os << boost::format("# %d x %d (0)\n") % n % n;
    
    for (int j=0; j<n; j++) {
      for (int i=0; i<n; i++)
	os << M(j, i) << " ";
      os << endl;
    }
    return(os);
  }

  template<class T>
  void writeTriangularAsFull(const string& fname, const Matrix<T, Triangular>& M, const string& meta) {
    ofstream ofs(fname.c_str());
    if (ofs == 0)
      throw(runtime_error("Cannot open " + fname + " for writing."));
    writeTriangularAsFull(ofs, M, meta);
  }


  // --- Reading ---

  template<class T, class P>
  istream& readAsciiMatrix(istream& is, Matrix<T,P>& M) {
    char inbuf[512];
    int m = 0, n = 0;

    // First, search for the marker...
    while (is.getline(inbuf, sizeof(inbuf))) {
      int i = sscanf(inbuf, "# %d x %d", &m, &n);
      if (i == 2)
	break;
    }

    if (m == 0 && n == 0)
      return(is);


    T datum;
    Matrix<T,P> R(m, n);
    for (int j=0;j < m; j++) {
      for (int i=0; i<n; i++) {
	if (!(is >> datum)) {
	  stringstream s;
	  s << "Invalid conversion on matrix read at (" << j << "," << i << ")";
	  throw(runtime_error(s.str()));
	}
	R(j, i) = datum;
      }
    }

    M = R;
    return(is);
  }

  template<class T, class P>
  void readAsciiMatrix(const string& fname, Matrix<T,P>& M) {
    ifstream ifs(fname.c_str());
    if (!ifs)
      throw(runtime_error("Cannot open " + fname + " for reading."));

    readAsciiMatrix(ifs, M);
  }

  // --- Specializations for reading ---

  template<class T>
  istream& readAsciiMatrix(istream& is, Matrix<T, Triangular>& M) {
    char inbuf[512];
    int m = 0;

    // First, search for the marker...
    while (is.getline(inbuf, sizeof(inbuf))) {
      char buf[20];
      int i = sscanf(inbuf, "# %d %10s", &m, buf);
      if (i != 2)
	continue;
      if (strncmp(buf, "TRIANGULAR", 10) == 0)
	break;
      throw(runtime_error("Magic matrix line found, but the matrix appears not to be triangular."));
    }

    if (m == 0)
      return(is);
    int n = m;

    Matrix<T, Triangular> R(m, n);
    long s = R.size();
    T datum;

    for (long i=0; i<s; i++) {
      if (!(is >> datum)) {
	stringstream s;
	s << "Invalid conversion on matrix read at [" << i << "]";
	throw(runtime_error(s.str()));
      }
      R[i] = datum;
    }

    M = R;
    return(is);
  }


  template<class T>
  void readAsciiMatrix(const string& fname, Matrix<T, Triangular>& M) {
    ifstream ifs(fname.c_str());
    if (!ifs)
      throw(runtime_error("Unable to open " + fname + " for reading."));

    readAsciiMatrix(ifs, M);
  }

}

#endif
