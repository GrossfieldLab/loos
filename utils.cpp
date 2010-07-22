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
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
#include <unistd.h>
#include <pwd.h>
#include <glob.h>

#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>

#include <boost/algorithm/string.hpp>

#include <Selectors.hpp>
#include <Parser.hpp>
#include <utils.hpp>



namespace loos {

  namespace {
    // Should this throw rather than exit?
    // Note that glob() expects this to return an int, but I think
    // it's better to abort with an error...
    int globError(const char* p, int errnum) {
      char *msg = strerror(errnum);
      std::cerr << msg << " - '" << p << "'" << std::endl;
      exit(-errnum);
    }
  };

  std::vector<std::string> globFilenames(const std::string& pattern) {
    glob_t buf;
    buf.gl_offs = 0;

    int i = glob(pattern.c_str(), 0, &globError, &buf);
    if (i != 0) {
      std::cerr << "Error code " << i << " from globbing '" << pattern << "'" << std::endl;
      exit(-255);
    }

    std::vector<std::string> names;
    for (uint j = 0; j<buf.gl_pathc; ++j)
      names.push_back(std::string( (buf.gl_pathv)[i] ));

    globfree(&buf);

    return(names);
  }


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

#if defined(REVISION)
    invoke += " [" + std::string(REVISION) + "]";
#endif

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
          throw(ParseError("Unable to parse " + s));
        if (!(is >> c.y()))
          throw(ParseError("Unable to parse " + s));
        if (!(is >> c.z()))
          throw(ParseError("Unable to parse " + s));

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
    catch(ParseError e) {
      throw(ParseError("Error in parsing '" + selection + "' ... " + e.what()));
    }

    KernelSelector selector(parser.kernel());
    AtomicGroup subset = source.select(selector);

    if (subset.size() == 0)
      throw(NullResult("No atoms were selected using '" + selection + "'"));

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


  template<>
  std::string parseStringAs<std::string>(const std::string& source, const uint pos, const uint nelem) {
    std::string val;

    uint n = !nelem ? source.size() - pos : nelem;
    if (pos + n > source.size())
      return(val);

    for (uint i=pos; i<pos+n; ++i)
      if (source[i] != ' ')
        val += source[i];

    return(val);
  }

  template<>
  std::string fixedSizeFormat(const std::string& s, const uint n) {
    uint m = s.size();
    if (m > n)
      return(s.substr(m-n, n));
    return(s);
  }


  namespace {
    int pow10[6]={1,10,100,1000,10000,100000};
    int pow36[6]={1,36,1296,46656,1679616,60466176};
  };


  
  int parseStringAsHybrid36(const std::string& source, const uint pos, const uint nelem) {
    uint n = !nelem ? source.size() - pos : nelem;
    if (pos + n > source.size())
      return(0);

    if (n > 6)
      throw(std::logic_error("Requested size exceeds max"));
    
    std::string s(source.substr(pos, n));
    bool negative(false);
    std::string::iterator si(s.begin());
    n = s.size();

    if (*si == '-') {
      negative = true;
      ++si;
      --n;
    }

    // Skip leading whitespace
    for (;si != s.end() && *si == ' '; ++si, --n) ;

    int offset = 0;   // This adjusts the range of the result
    char cbase = 'a'; // Which set or characters (upper or lower) for the alpha-part
    int ibase = 10;   // Number-base (i.e. 10 or 36)

    // Decide which chunk we're in...
    if (*si >= 'a') {
      offset = pow10[n] + 16*pow36[n-1];
      cbase = 'a';
      ibase = 36;
    } else if (*si >= 'A') {
      offset = pow10[n] - 10*pow36[n-1];
      cbase = 'A';
      ibase = 36;
    }

    int result = 0;
    while (si != s.end()) {
      int c = (*si >= cbase) ? *si-cbase+10 : *si-'0';
      result = result * ibase + c;
      ++si;
    }

    result += offset;

    if (negative)
      return(-result);

    return(result);
  }



  // Note: this currently will overflow if sufficiently negative
  // to overflow the base-10 part...
  std::string hybrid36AsString(int d, uint n) {

    if (n > 6)
      throw(std::logic_error("Requested size exceeds max"));

    int n10 = pow10[n];
    int n36 = pow36[n-1];
    int cuta = n10 + n36 * 26;  // Cutoff between upper and lower
                                // representations (i.e. A000 vs a000)
    bool negative(false);

    if (d < 0) {
      negative = true;
      d = -d;
    }

    if (d >= n10 + 52 * n36)
      throw(std::runtime_error("Number out of range"));

    unsigned char coffset = '0';   // Digits offset for output
    int ibase = 10;                // Numeric base (i.e. 10 or 36)

    if (d >= cuta) {
      coffset = 'a' - 10;
      ibase = 36;
      d -= cuta;
      d += 10*n36;
    } else if (d >= n10) {
      coffset = 'A' - 10;
      d -= n10;
      d += 10*n36;
      ibase = 36;
    }

    std::string result;
    while (d > 0) {
      unsigned char digit = d % ibase;
      digit += (digit > 9) ? coffset : '0';

      result.push_back(digit);
      d /= ibase;
    }

    if (negative)
      result.push_back('-');

    // right-justify...should we be left instead??
    for (uint i=result.size(); i < n; ++i)
      result.push_back(' ');

    std::reverse(result.begin(), result.end());
    return(result);
  }



}
