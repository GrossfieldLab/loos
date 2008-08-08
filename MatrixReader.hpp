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



#if !defined(MATRIXREADER_HPP)
#define MATRIXREADER_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <Fmt.hpp>
#include <string>
#include <cctype>

#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/string.hpp>



using namespace std;



//! Base class for reading matrices
template<class T>
class MatrixReader {
public:
  typedef unsigned int uint;
  typedef boost::tuple<T*, unsigned int, unsigned int> Result;

  MatrixReader() : ins(&cin) { }
  explicit MatrixReader(istream* i) : ins(i) { }
  virtual ~MatrixReader() { }


  virtual Result basic_read(istream* input) =0;

  virtual Result read(istream* input) {
    return(basic_read(input));
  }

  virtual Result read(const string& s) {
    ifstream ifs(s.c_str());
    if (!ifs)
      throw(runtime_error("Unable to open " + s));
    return(MatrixReader<T>::read(&ifs));
  }


  virtual Result read(const char *s) {
    ifstream ifs(s);
    if (!ifs)
      throw(runtime_error("Unable to open " + string(s)));
    return(MatrixReader<T>::read(&ifs));
  }

protected:
  istream* ins;
};


//! Class for reading raw ascii matrices
/** This reader will skip over any lines that do not begin with a
 * digit at the start of the file.  This is the associated metadata
 * (which is currently unused).  Reading the matrix continues until a
 * line without a digit at the start.
 *
 * Note that although the matrix is expected to be written in
 * row-major order, it is actually placed in memory as col-major order
 * since we're assuming it'll be fed to lapack/blas...
 */

template<class T>
class RawAsciiReader : public MatrixReader<T> {
public:
  static const uint inbufsiz = 256536;

  typename MatrixReader<T>::Result basic_read(istream* input);

private:
  boost::tuple<uint, uint> scanSize(istream* input);
};



template<class T> boost::tuple<uint,uint> RawAsciiReader<T>::scanSize(istream* input) {

  uint m, n;
  char inbuf[inbufsiz];

  // Record the start of the raw matrix data...
  streampos curpos = input->tellg();
  int i;
  while (input->getline(inbuf, inbufsiz)) {
    stringstream sin(inbuf);
    if (sin >> i)
      break;
    curpos = input->tellg();
  }

  // Got the first line, count the number of columns...
  vector<string> strings;
  string line(inbuf);
  boost::trim(line);
  boost::split(strings, line, boost::is_any_of(" \t"), boost::token_compress_on);
  n = strings.size();
  if (n == 0)
    throw(runtime_error("Could not find any columns in the matrix!"));
  // Now count the rows...
  m=1;
  while (input->getline(inbuf, inbufsiz)) {
    stringstream sin(inbuf);
    if (!(sin >> i)) {
      cerr << "Error while reading ASCII matrix at " << inbuf[0] << inbuf[1] << inbuf[2] << endl;
      break;
    }
    ++m;
  }

  // Restore filepos
  input->clear();
  input->seekg(curpos);
  boost::tuple<uint, int> result(m, n);

  return(result);
}


template<class T> typename MatrixReader<T>::Result RawAsciiReader<T>::basic_read(istream* input) {

  boost::tuple<uint,uint> size = scanSize(input);
  uint m = boost::get<0>(size);
  uint n = boost::get<1>(size);

  // Should we impose a max matrix size here???
  T* data = new T[m*n];
  T datum;

  uint i, j;
  for (j=0; j<m; j++)
    for (i=0; i<n; i++) {
      if (!(*input >> datum)) {
	stringstream s;
	s << "Invalid conversion on matrix read at (" << j << "," << i << ") [" << datum << "]";
	throw(runtime_error(s.str()));
      }
      data[i*m+j] = datum;
    }

  typename MatrixReader<T>::Result result(data, m, n);
  return(result);
}



#endif
