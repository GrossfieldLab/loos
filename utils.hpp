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


#include "Coord.hpp"

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

#endif
