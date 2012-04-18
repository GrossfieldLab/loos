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




#if !defined(LOOS_UTILS_HPP)
#define LOOS_UTILS_HPP

#include <iostream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <exception>
#include <stdexcept>


#include <boost/algorithm/string.hpp>
#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include <ctime>


#include <loos_defs.hpp>
#include <exceptions.hpp>
#include <Coord.hpp>
#include <pdb_remarks.hpp>



//! Namespace for most things not already encapsulated within a class.
namespace loos {

  //! Pull off the file name extension (if present)
  std::string findBaseName(const std::string&);

  //! Get the next line of input, skipping blanks and stripping comments
  std::string getNextLine(std::istream& is, int* lineno);

  //! Read a list of numbers from a stream
  template<typename T>
  std::vector<T> readVector(std::istream& is) {
    std::vector<T> data;
    for (;;) {
      std::string s = getNextLine(is, 0);
      if (s.length() == 0)
        break;

      std::istringstream iss(s);
      T datum;
      iss >> datum;
      data.push_back(datum);
    }

    return(data);
  }

  template<typename T>
  std::vector< std::vector<T> > readTable(std::istream& is) {
    std::vector< std::vector<T> > table;
    for (;;) {
      std::string s = getNextLine(is, 0);
      if (s.empty())
        break;

      std::istringstream iss(s);
      T datum;
      std::vector<T> row;
      while (iss >> datum)
        row.push_back(datum);
      table.push_back(row);
    }
    return(table);
  }

  //! Create an invocation header
  /**
   *This is a string that can be embedded in output that records the
   *invoking user, command-line, and a timestamp.
   */
  std::string invocationHeader(int, char *[]);





  //! Parse an Octave/Matlab-style range
  /**
   * The Octave format is one of the following:
   * - value
   * - start:stop
   * - start:step:stop
   *
   * The range is inclusive for both ends.  Internally, a vector of
   * T's is created.  There is no checking to make sure that the
   * vector doesn't completely fill up memory.
   *
   * As with matlab/octave, to count down, the step must be negative
   * and start > stop.
   *
   */
  template<typename T>
  std::vector<T> parseRange(const std::string& text) {
    T a;
    T b;
    T c;
    char sep;
    std::vector<T> indices;
    std::istringstream is(text);
    
    is >> a;
    if (is.eof()) {
      indices.push_back(a);
      return(indices);
    }

    is >> sep;
    bool is_negative = false;
    if (is.peek() == '-') {
      is_negative = true;
      is.get();
    }
    is >> b;
    if (is.fail() || sep != ':')
      throw(ParseError("Could not parse range " + text));
    
    if (is.eof()) {
      c = 1;
      if (is_negative)
        b = -b;
      
      is_negative = a > b;
        
    } else {
      c = b;
        is >> sep >> b;
        if (is.fail() || sep != ':')
          throw(ParseError("Could not parse range " + text));
    }

    // Some input validation...
    if (a > b && !is_negative)
      throw(ParseError("You must use a negative step to count down: " + text));
    if (a < b && is_negative)
      throw(ParseError("You must use a postive step to count up: " + text));
    if (c == 0)
      throw(ParseError("Thou shalt only use non-zero step sizes: " + text));

    if (is_negative) {

      T i;
      for (i=a; i > b; i -= c)
        indices.push_back(i);
      if (a < -a && b == 0)     // If unsigned type, cannot use >= 0 as test
        indices.push_back(0);
      else if (i >= b)
        indices.push_back(i);

    } else
      for (T i=a; i <= b; i += c)
        indices.push_back(i);

    return(indices);
  }




  //! Parses a comma-separated list of Octave-style ranges
  /**
   * This function breaks apart a string at the commas and passes each
   * substring to parseRange.  The union of all of the vectors
   * returned by parseRange is then sorted in ascending order and
   * duplicate values are removed.  This vector is then returned to
   * the caller.
   */
  template<typename T>
  std::vector<T> parseRangeList(const std::string& text) {
    std::vector<std::string> terms;
    std::set<T> indices;
    std::insert_iterator< std::set<T> > ii(indices, indices.begin());

    boost::split(terms, text, boost::is_any_of(","), boost::token_compress_on);
    std::vector<std::string>::const_iterator ci;
    for (ci = terms.begin(); ci != terms.end(); ci++) {
      if (ci->empty())
        continue;
      std::vector<T> result = parseRange<T>(*ci);
      std::copy(result.begin(), result.end(), ii);
    }
    std::vector<T> results(indices.size());
    std::copy(indices.begin(), indices.end(), results.begin());
    return(results);
  }


  //! Parses a list of Octave-style range specifiers (for compatability)
  std::vector<int> parseRangeList(const std::string&);

  //! Parses a list of Octave-style ranges taken from a vector of strings
  template<typename T>
  std::vector<T> parseRangeList(const std::vector<std::string>& ranges) {
    std::ostringstream os;
    std::copy(ranges.begin(), ranges.end(), std::ostream_iterator<std::string>(os, ","));
      return(parseRangeList<T>(os.str()));
  }

  //! Applies a string-based selection to an atomic group...
  AtomicGroup selectAtoms(const AtomicGroup&, const std::string);


  //! Returns a byte-swapped copy of an arbitrary type
  /** 
   * Only valid for simple types (i.e. int, float, double)
   */
  template<typename T>
  T swab(const T& datum) {
    uint size = sizeof(T);
    const unsigned char* p = reinterpret_cast<const unsigned char*>(&datum);
    T swabbed;
    unsigned char* q = reinterpret_cast<unsigned char*>(&swabbed);
    
    uint i, j;
    for (i=0, j=size-1; i<size; ++i, --j)
      q[i] = p[j];
    
    return(swabbed);
  }

  std::string timeAsString(const double t);

  //! Extracts a field from a string
  template<typename T>
  T parseStringAs(const std::string& source, const uint pos =0, const uint nelem =0) {
    T val(0);

    uint n = !nelem ? source.size() - pos : nelem;
    if (pos + n > source.size())
      return(val);
    
    std::string element(source.substr(pos, n));
    std::istringstream iss(element);
    if (!(iss >> val)) {
      std::stringstream msg;
      msg << "PARSE ERROR\n" << source << std::endl;
      for (uint i=0; i<pos; ++i)
        msg << ' ';
      msg << '^';
      if (n > 1)
        for (uint i=1; i<n; ++i)
          msg << '^';
      msg << std::endl;
      throw(ParseError(msg.str()));
    }

    return(val);
  }

  template<> std::string parseStringAs<std::string>(const std::string& source, const uint pos, const uint nelem);

  template<typename T>
  std::string fixedSizeFormat(const T t, const uint n) {
    std::stringstream ss;
    ss << t;
    std::string s(ss.str());
    uint m = s.size();
    if (m > n)
      return(s.substr(m-n, n));
    return(s);
  }

  template<> std::string fixedSizeFormat(const std::string& s, const uint n);

  //! Convert a hybrid-36 encoded string into an int
  int parseStringAsHybrid36(const std::string& source, const uint pos =0, const uint nelem =0);

  //! Convert an int into a hybrid-36 encoded string
  std::string hybrid36AsString(int value, uint fieldsize);

  // The following are for support of boost::program_options

  //! Convert something that can iterate into a string...
  template<typename T> std::string vToString(const T& x) {
    std::ostringstream oss;

    for (typename T::const_iterator i = x.begin(); i != x.end(); ++i)
      oss << *i << ((i == x.end() - 1) ? "" : ",");

    return(oss.str());
  }


  //! Removes internal newlines from string
  std::string sanitizeString(const std::string& s);

  //! Converts a vector of strings into a standard log format
  std::string stringsAsComments(const std::vector<std::string>& v);


  //! Converts a vector of strings into a single string with newlines
  std::string stringsAsString(const std::vector<std::string>& v);



  //! Convert a vector of type T to a string-list with commas
  template<typename T> std::string vectorAsStringWithCommas(const std::vector<T>& v) {
    std::ostringstream oss;
    for (typename std::vector<T>::const_iterator i = v.begin(); i != v.end(); ++i) {
      oss << *i;
      if (i != v.end() - 1)
        oss << ",";
    }
    return(oss.str());
  }



};

#endif
