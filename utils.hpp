/*
  utils.hpp
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Suite-wide utilities...
*/


#if !defined(UTILS_HPP)
#define UTILS_HPP

#include <iostream>
#include <vector>


#include <boost/random.hpp>
#include <ctime>


#include <loos.hpp>
#include <Coord.hpp>
#include <pdb_remarks.hpp>

using namespace std;
using namespace __gnu_cxx;

//! Get the next line of input, skipping blanks and stripping comments
string getNextLine(istream&, int*);

//! Read a list of integers from a stream
vector<int> readIndexMap(istream&);

//! Create an invocation header
/*!
  This is a string that can be embedded in output that records the
  invoking user, command-line, and a timestamp.
*/
string invocationHeader(int, char *[]);

//! Extract the Alan-style box-size from a PDB Remarks block.
/** Returns a GCoord(99999.99, 99999.99, 99999.99) if there is no box
 *  info found in the remarks block.
 */
GCoord boxFromRemarks(const Remarks&);

//! Checks to see if a Remarks block has an Alan-style box size in it.
bool remarksHasBox(const Remarks&);

// The following are in LOOS namespace because they are either
// collisions waiting to or are too esoteric to warrant going into std

namespace loos {


  typedef boost::mt19937 base_generator_type;

  //! Suite-wide random number generator singleton
  base_generator_type& rng_singleton(void);

  //! Randomly seeds the RNG
  /**Currently uses time(3) to seed the RNG obtained from the singleton...
   */
  void randomSeedRNG(void);

};

#endif
