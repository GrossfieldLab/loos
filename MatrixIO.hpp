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
#include <Fmt.hpp>
#include <string>


using namespace std;

class MatrixWriter {

protected:
  typedef unsigned int uint;

  class Wrapper {
  public:
    virtual double operator()(const uint i) const =0;
    virtual ~Wrapper() { }
  };

  class FloatWrapper : public Wrapper {
  public:
    FloatWrapper(const float *p) : data(p) { }
    double operator()(const uint i) const { return(data[i]); }
  private:
    const float *data;
  };

  class DoubleWrapper : public Wrapper {
  public:
    DoubleWrapper(const double *p) : data(p) { }
    double operator()(const uint i) const { return(data[i]); }
  private:
    const double *data;
  };



public:

  MatrixWriter() : prefixname(""), ofs(&cout) { }
  MatrixWriter(const string& prefix) : prefixname(prefix), ofs(0) { }
  MatrixWriter(ostream* ofsp) : prefixname(""), ofs(ofsp) { }
  virtual ~MatrixWriter() { }
  

  string prefix(void) const { return(prefixname); }
  void prefix(const string& s) { prefixname = s; }
  ostream* outputStream(void) const { return(ofs); }
  void outputStream(ostream* ofsp) { ofs = ofsp; }
  string metadata(void) const { return(meta_data); }
  void metadata(const string& s) { meta_data = s; }
  

  void basic_write(const Wrapper& data, const string& tag, const uint m, const uint n, const bool trans, const uint maxcol, const uint maxrow);

  void write(const float* p, const string& tag, const uint m, const uint n, const bool trans = false, const uint maxcol=0, const uint maxrow=0) {
    basic_write(FloatWrapper(p), tag, m, n, trans, maxcol, maxrow);
  }

  void write(const double* p, const string& tag, const uint m, const uint n, const bool trans = false, const uint maxcol=0, const uint maxrow=0) {
    basic_write(DoubleWrapper(p), tag, m, n, trans, maxcol, maxrow);
  }

  virtual void OutputPreamble(ostream *po, const string& tag, const uint m, const uint n, const bool trans) =0;
  virtual void OutputDatum(ostream *po, const double d) =0;
  virtual void OutputEOL(ostream *po) =0;
  virtual void OutputCoda(ostream *po) =0;
  virtual string constructFilename(const string& tag) =0;

protected:
  string prefixname;
  string meta_data;
  ostream *ofs;
};



class RawAsciiWriter : public MatrixWriter {
public:
  RawAsciiWriter() : MatrixWriter() { }
  RawAsciiWriter(const string& s) : MatrixWriter(s) { }
  RawAsciiWriter(ostream* o) : MatrixWriter(o) { }

  void OutputPreamble(ostream *po, const string& tag, const uint m, const uint n, const bool trans) {
    if (meta_data != "")
      *po << "# " << meta_data << endl;
    *po << "# " << m << " " << n << " " << trans << " \"" << tag << "\"" << endl;
  }

  void OutputDatum(ostream *po, const double d) {
    *po << d << " ";
  }

  void OutputEOL(ostream *po) {
    *po << endl;
  }
  
  void OutputCoda(ostream *po) { }

  string constructFilename(const string& tag) { return(string(prefixname + tag + ".asc")); }
};



class OctaveAsciiWriter : public MatrixWriter {
public:
  OctaveAsciiWriter() : MatrixWriter() { }
  OctaveAsciiWriter(const string& s) : MatrixWriter(s) { }
  OctaveAsciiWriter(ostream* o) : MatrixWriter(o) { }

  void OutputPreamble(ostream *po, const string& tag, const uint m, const uint n, const bool trans) {
    if (meta_data != "")
      *po << "% " << meta_data << endl;
    *po << tag << " = [\n";
  }

  void OutputDatum(ostream *po, const double d) {
    *po << d << " ";
  }

  void OutputEOL(ostream *po) {
    *po << " ;\n";
  }

  void OutputCoda(ostream *po) {
    *po << "];\n";
  } 

  string constructFilename(const string& tag) { return(string(prefixname + tag + ".m")); }

};



#endif
