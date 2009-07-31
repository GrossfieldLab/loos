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




#if !defined(UTILS_HPP)
#define UTILS_HPP

#include <iostream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <stdexcept>


#include <boost/random.hpp>
#include <boost/algorithm/string.hpp>

#include <ctime>


#include <loos_defs.hpp>
#include <Coord.hpp>
#include <pdb_remarks.hpp>


//! Namespace for most things not already encapsulated within a class.
namespace loos {

  //! Pull off the file name extension (if present)
  std::string findBaseName(const std::string&);

  //! Get the next line of input, skipping blanks and stripping comments
  std::string getNextLine(std::istream&, int*);

  //! Read a list of integers from a stream
  std::vector<int> readIndexMap(std::istream&);

  //! Create an invocation header
  /**
   *This is a string that can be embedded in output that records the
   *invoking user, command-line, and a timestamp.
   */
  std::string invocationHeader(int, char *[]);

  //! Extract the Alan-style box-size from a PDB Remarks block.
  /** Returns a GCoord(99999.99, 99999.99, 99999.99) if there is no box
   *  info found in the remarks block.
   */
  GCoord boxFromRemarks(const Remarks&);

  //! Checks to see if a Remarks block has an Alan-style box size in it.
  bool remarksHasBox(const Remarks&);

  // The following are in LOOS namespace because they are either
  // collisions waiting to or are too esoteric to warrant going into std

  typedef boost::mt19937 base_generator_type;

  //! Suite-wide random number generator singleton
  base_generator_type& rng_singleton(void);

  //! Randomly seeds the RNG
  /**Currently uses time(3) to seed the RNG obtained from the singleton...
   */
  void randomSeedRNG(void);

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
   * Care should also be exercised when using unsigned types and
   * reversed ranges, i.e. parseRange<unsigned int>("5:0")
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

    is >> sep >> b;
    if (is.fail() || sep != ':')
      throw(std::runtime_error("Could not parse range " + text));
    
    if (is.eof())
      c = 1;
    else {
      c = b;
        is >> sep >> b;
        if (is.fail() || sep != ':')
          throw(std::runtime_error("Could not parse range " + text));
    }

    // Catch a potentially sticky situation...
    if (a < -a && b == 0)
      throw(std::logic_error("Using an unsigned type with a reverse-range ending at zero."));

    if (b >= a)
      for (T i=a; i <= b; i += c)
        indices.push_back(i);
    else
      for (T i=a; i >= b; i -= c)
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
    T val;

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
      throw(std::runtime_error(msg.str()));
    }

    return(val);
  }

  template<> std::string parseStringAs<std::string>(const std::string& source, const uint pos, const uint nelem);

};

#endif
