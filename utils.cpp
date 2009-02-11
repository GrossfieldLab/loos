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
#include <iomanip>

#include <boost/algorithm/string.hpp>

#include <Selectors.hpp>
#include <Parser.hpp>
#include <utils.hpp>


namespace loos {

  
  std::string findBaseName(const std::string& s) {
    std::string result;

    int n = s.find('.');
    result = (n <= 0) ? s : s.substr(0, n);
    
    return(result);
  }


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



  base_generator_type& rng_singleton(void) {
    static base_generator_type rng;

    return(rng);
  }


  void randomSeedRNG(void) {
    base_generator_type& rng = rng_singleton();

    rng.seed(static_cast<unsigned int>(time(0)));
  }


  std::vector<int> parseRangeList(const std::string& text) {
    return(parseRangeList<int>(text));
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
  AtomicGroup selectAtoms(const AtomicGroup& source, const std::string selection) {
  
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

  std::string timeAsString(const double t) {
    if (t < 90.0) {
      std::stringstream s;
      s << std::fixed << std::setprecision(3) << t << "s";
      return(s.str());
    }
  
    double mins = floor(t / 60.0);
    double secs = t - mins * 60.0;
    if (mins < 90.0) {
      std::stringstream s;
      s << std::fixed << std::setprecision(0) << mins << "m" << std::setprecision(3) << secs << "s";
      return(s.str());
    }
  
    double hrs = floor(mins / 60.0);
    mins -= hrs * 60.0;
    std::stringstream s;
    s << std::fixed << std::setprecision(0) << hrs << "h" << mins << "m" << std::setprecision(3) << secs << "s";
    return(s.str());
  }


}
