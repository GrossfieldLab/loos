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



#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pwd.h>

#include <algorithm>
#include <string>
#include <sstream>

#include <boost/algorithm/string.hpp>

#include <Selectors.hpp>
#include <Parser.hpp>
#include <utils.hpp>


std::string getNextLine(std::istream& is, int *lineno = 0) {
  std::string s;
  std::string::size_type i;

  for (;;) {
    if (getline(is, s, '\n').eof())
      break;

    if (lineno != 0)
      *lineno += 1;

    // Strip off comments
    i = s.find('#');
    if (i != std::string::npos)
      s.erase(i, s.length() - i);

    // Remove leading whitespace
    i = s.find_first_not_of(" \t");
    if (i > 0)
      s.erase(0, i);

    // Is std::string non-empty?
    if (s.length() != 0)
      break;
  }

  return(s);
}



std::vector<int> readIndexMap(std::istream& is) {
  std::vector<int> indices;
  int i;
  std::string s;

  for (;;) {
    s = getNextLine(is);
    if (s.length() == 0)
      break;

    std::istringstream iss(s);
    iss >> i;
    indices.push_back(i);
  }

  return(indices);
}



std::string invocationHeader(int argc, char *argv[]) {
  std::string invoke, user;
  int i;
  
  time_t t = time(0);
  std::string timestamp(asctime(localtime(&t)));
  timestamp.erase(timestamp.length() - 1, 1);

  struct passwd* pwd = getpwuid(getuid());
  if (pwd == 0)
    user = "UNKNOWN USER";
  else
    user = pwd->pw_name;

  invoke = std::string(argv[0]) + " ";
  std::string sep(" ");
  for (i=1; i<argc; i++) {
    if (i == argc-1)
      sep = "";
    invoke += "'" + std::string(argv[i]) + "'" + sep;
  }

  invoke += " - " + user + " (" + timestamp + ")";

  // Since some args my be brought in from a file via the shell
  // back-tick operator, we process embedded returns...
  boost::replace_all(invoke, "\n", "\\n");

  return(invoke);
}



GCoord boxFromRemarks(const Remarks& r) {
  int n = r.size();
  int i;

  GCoord c(99999.99, 99999.99, 99999.99);

  for (i=0; i<n; i++) {
    std::string s = r[i];
    if (s.substr(0, 6) == " XTAL ") {
      std::stringstream is(s.substr(5));
      if (!(is >> c.x()))
        throw(std::runtime_error("Unable to parse " + s));
      if (!(is >> c.y()))
        throw(std::runtime_error("Unable to parse " + s));
      if (!(is >> c.z()))
        throw(std::runtime_error("Unable to parse " + s));

      break;
    }
  }

  return(c);
}



bool remarksHasBox(const Remarks& r) {
  int n = r.size();
  for (int i = 0; i<n; i++) {
    std::string s = r[i];
    if (s.substr(0, 6) == " XTAL ")
      return(true);
  }
  return(false);
}



loos::base_generator_type& loos::rng_singleton(void) {
  static base_generator_type rng;

  return(rng);
}


void loos::randomSeedRNG(void) {
  base_generator_type& rng = rng_singleton();

  rng.seed(static_cast<unsigned int>(time(0)));
}

/** This routine breaks the input string into chunks delimited by
 * commas.  Each chunk is a range specified in Octave format, i.e.
 * - start:stop
 * - start:step:stop
 *
 * The range is inclusive of both ends.  Additionally, a single index
 * can be specified.
 *
 * Internally, this routine creates a vector of ints that represent
 * the specified indices.  There is no bounds checking...  Duplicate
 * indices are filtered and the returned vector is sorted.
 */ 
std::vector<int> loos::parseRangeList(const std::string& text) {
  std::vector<std::string> terms;
  std::vector<int> indices;

  boost::split(terms, text, boost::is_any_of(","), boost::token_compress_on);
  std::vector<std::string>::const_iterator ci;
  for (ci = terms.begin(); ci != terms.end(); ci++) {
    int a, b, c;
    int i;
    i = sscanf(ci->c_str(), "%d:%d:%d", &a, &b, &c);
    if (i == 2) {
      c = b;
      b = 1;
    } else if (i == 1) {
      c = a;
      b = 1;
    } else if (i != 3)
      throw(std::runtime_error("Cannot parse range list item " + *ci));

    if (c < a) {
      if (b > 0)
        throw(std::runtime_error("Invalid range spec " + *ci));
      int x = c;
      c = a;
      a = x;
      b = -b;
    } else if (b <= 0)
      throw(std::runtime_error("Invalid range spec " + *ci));

    for (int i=a; i<=c; i += b)
      indices.push_back(i);
  }
  sort(indices.begin(), indices.end());
  std::vector<int> results;
  std::vector<int>::const_iterator cvi;
  int last = indices[0];
  results.push_back(last);

  for (cvi = indices.begin()+1; cvi != indices.end(); cvi++)
    if (*cvi != last) {
      last = *cvi;
      results.push_back(last);
    }

  return(results);
}


/** This routine parses the passed string, turning it into a selector
 *  and applies it to \a source.  If there is an exception in the
 *  parsing, this is repackaged into a more sensible error message
 *  (including the string that generated the error).  No other
 *  exceptions are caught.
 *
 *  We're also assuming that you're <EM>always</EM> wanting to select
 *  some atoms, so lack of selection constitutes an error and an
 *  exception is thrown.  Note that in both the case of a parse error
 *  and null-selection, a runtime_error exception is thrown so the
 *  catcher cannot disambiguate between the two.
*/
AtomicGroup loos::selectAtoms(const AtomicGroup& source, const std::string selection) {
  
  Parser parser;

  try {
    parser.parse(selection);
  }
  catch(std::runtime_error e) {
    throw(std::runtime_error("Error in parsing '" + selection + "' ... " + e.what()));
  }

  KernelSelector selector(parser.kernel());
  AtomicGroup subset = source.select(selector);

  if (subset.size() == 0)
    throw(std::runtime_error("No atoms were selected using '" + selection + "'"));

  return(subset);
}
