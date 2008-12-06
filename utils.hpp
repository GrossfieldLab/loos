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
#include <string>
#include <vector>
#include <stdexcept>


#include <boost/random.hpp>
#include <boost/algorithm/string.hpp>
#include <ctime>


#include <loos_defs.hpp>
#include <Coord.hpp>
#include <pdb_remarks.hpp>


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

//! Namespace for segregating esoteric functions or functions with common names
namespace loos {


  typedef boost::mt19937 base_generator_type;

  //! Suite-wide random number generator singleton
  base_generator_type& rng_singleton(void);

  //! Randomly seeds the RNG
  /**Currently uses time(3) to seed the RNG obtained from the singleton...
   */
  void randomSeedRNG(void);

  //! Parses a list of Octave-style range specifiers
  std::vector<int> parseRangeList(const std::string& text);

  //! Applies a string-based selection to an atomic group...
  AtomicGroup selectAtoms(const AtomicGroup&, const std::string);

};

#endif
