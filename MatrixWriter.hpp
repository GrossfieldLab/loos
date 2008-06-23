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
#include <Fmt.hpp>
#include <string>


using namespace std;


//! Class for handling writing of matrices in different formats.
template<class T>
class MatrixWriter {
public:
  typedef unsigned int uint;

  //! Default constructor sends output to cout
  MatrixWriter() : prefixname(""), ofs(&cout) { }

  //! If passed a string, this forces output to open files named prefix + tag
  explicit MatrixWriter(const string& prefix) : prefixname(prefix), ofs(0) { }
  
  //! If passed a pointer to an ostream, then that is used instead...
  explicit MatrixWriter(ostream* ofsp) : prefixname(""), ofs(ofsp) { }
  virtual ~MatrixWriter() { }
  

  string prefix(void) const { return(prefixname); }
  void prefix(const string& s) { prefixname = s; }
  ostream* outputStream(void) const { return(ofs); }
  void outputStream(ostream* ofsp) { ofs = ofsp; }

  string metadata(void) const { return(meta_data); }
  //! Metadata is written out (optionally) depending on format.
  void metadata(const string& s) { meta_data = s; }
  

  //! This handles writing the matrix out.
  /** It is assumed that matrices are stored col-major (i.e. Fortran-style)
   *  Notable parameters:
   *  \arg \c tag String used to tag the matrix (i.e. name it)
   *  \arg \c trans The matrix is transposed on output
   *  \arg \c maxcol The maximum column to write
   *  \arg \c maxrow The maxmimum row to write
   */
  void write(const T* data, const string& tag, const uint m, const uint n, const bool trans = false, const uint maxcol = 0, const uint maxrow = 0);


  // These are overriden by subclasses to control the output format...

  //! Writes out any pre-matrix text as required by the format.
  virtual void OutputPreamble(ostream *po, const string& tag, const uint m, const uint n, const bool trans) =0;

  //! Writes out a single element of the matrix
  virtual void OutputDatum(ostream *po, const T d) =0;

  //! Ends a row (line) of data
  virtual void OutputEOL(ostream *po) =0;

  //! Anything that needs to come after the matrix data has been written...
  virtual void OutputCoda(ostream *po) =0;

  //! Constructs format-dependent output filename...
  virtual string constructFilename(const string& tag) =0;

protected:
  string prefixname;
  string meta_data;
  ostream *ofs;
};


//! Class for writing raw ascii matrices
/** Matrix properties (such as size and transpose flags) are written
 *  in the preamble
 */
template<class T>
class RawAsciiWriter : public MatrixWriter<T> {
public:
  RawAsciiWriter() : MatrixWriter<T>() { }
  explicit RawAsciiWriter(const string& s) : MatrixWriter<T>(s) { }
  explicit RawAsciiWriter(ostream* o) : MatrixWriter<T>(o) { }

  void OutputPreamble(ostream *po, const string& tag, const uint m, const uint n, const bool trans) {
    if (RawAsciiWriter<T>::meta_data != "")
      *po << "# " << RawAsciiWriter<T>::meta_data << endl;
    *po << "# " << m << " " << n << " " << trans << " \"" << tag << "\"" << endl;
  }

  void OutputDatum(ostream *po, const T d) {
    *po << d << " ";
  }

  void OutputEOL(ostream *po) {
    *po << endl;
  }
  
  void OutputCoda(ostream *po) { }

  string constructFilename(const string& tag) { return(string(RawAsciiWriter<T>::prefixname + tag + ".asc")); }
};


//! Class for writing ASCII Octave format (as in a .m script)
template<class T>
class OctaveAsciiWriter : public MatrixWriter<T> {
public:
  OctaveAsciiWriter() : MatrixWriter<T>() { }
  explicit OctaveAsciiWriter(const string& s) : MatrixWriter<T>(s) { }
  explicit ctaveAsciiWriter(ostream* o) : MatrixWriter<T>(o) { }

  void OutputPreamble(ostream *po, const string& tag, const uint m, const uint n, const bool trans) {
    if (OctaveAsciiWriter<T>::meta_data != "")
      *po << "% " << OctaveAsciiWriter<T>::meta_data << endl;
    *po << tag << " = [\n";
  }

  void OutputDatum(ostream *po, const T d) {
    *po << d << " ";
  }

  void OutputEOL(ostream *po) {
    *po << " ;\n";
  }

  void OutputCoda(ostream *po) {
    *po << "];\n";
  } 

  string constructFilename(const string& tag) { return(string(OctaveAsciiWriter<T>::prefixname + tag + ".m")); }

};





template<class T>
void MatrixWriter<T>::write(const T* data, const string& tag, const uint m, const uint n, const bool trans, const uint maxcol, const uint maxrow) {
  uint i, j, k, s = m*n;
  uint nn = (maxcol == 0 || maxcol > n) ? n : maxcol;
  uint mm = (maxrow == 0 || maxrow > n) ? m : maxrow;
  ofstream *ofsp = 0;
  
  // Determine if we need to open a file or not...
  ostream *po = ofs;
  if (!po) {
    string fname = constructFilename(tag);
    ofsp = new ofstream(fname.c_str());
    if (!ofsp)
      throw(runtime_error("Unable to open file " + fname));
    po = ofsp;
  }

  OutputPreamble(po, tag, m, n, trans);
  for (j=0; j<mm; j++) {
    for (i=0; i<nn; i++) {
      double d;

      if (trans)
	k = j*n+i;
      else
	k = i*m+j;
      assert(k < s);
      d = data[k];
      OutputDatum(po, d);
    }
    OutputEOL(po);
  }
  OutputCoda(po);
  
  delete ofsp;
}



#endif
