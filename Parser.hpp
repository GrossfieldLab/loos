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




#if !defined(PARSER_HPP)
#define PARSER_HPP

#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

#include <AtomicGroup.hpp>
#include <Kernel.hpp>
#include <ParserDriver.hpp>


//! Front-end to the Bison/Flex parser
/** Creates a compiled Kernel that can then be executed to select
 *  atoms.  The grammar accepted is relatively simple and patterned
 *  after C/PERL expressions.  Relational operators are allowed, as
 *  are basic logical operators (and, or, and not).  Pre-defined
 *  keywords are: name, id, resname, resid, segid.  These evaluate to
 *  the current atom's appropriate property.  Case IS significant both
 *  for keywords and for strings.  Integer numbers are allowed.
 *  Strings are delimited by either single quotes (') or double quotes
 *  (").  Also, strings inequalities are handled lexically as C++
 *  strings are, i.e. string1 > string2 is true is it would be in
 *  regular C++ code.
 *
 *  Regular expressions (in PERL format) are supported via Boost.  The
 *  regular expression matching operator, '=~' is slightly special in
 *  that it will only permit you to match a keyword that would
 *  evaluate to a string.  In other words, you may match against a
 *  name, resname, and segname (segname is an alias for segid), but
 *  NOT an id nor a resid. 
 *
 *  String equality matches the entire string.  If you want to match a
 *  subset, you should use the '=~' operator.  In other words,
    \verbatim
    "CA" == "C"   -> false
    "C"  == "C"   -> true
    "CA" =~ "C"   -> true
    \endverbatim
 *
 *  Finally, the standard precedence and associativity that apply in
 *  C++ apply here.  Expressions are evaluated left to right and
 *  parenthesis may be used to alter precedence/evaluation order.
 *  Unlike C/C++ however, the logical operators do not short-circuit.
 *
 *  If there is a syntax error in the selection string, then a
 *  runtime_error() is thrown.
 *
 *  Examples:
 *  \code
 *  string selection_string = "resid >= 10 && resid <= 100 && name == 'CA'";
 *  Parser parsed(selection_string);
 *  KernelSelector parsed_selector(parsed.kernel());
 *  AtomicGroup parsed_selection = molecule.select(parsed_selector)
 *  \endcode
 *
 *  Parser objects are intended to be a parse-once object.  If you
 *  want to parse multiple selection strings, then you should
 *  instantiate a Parser object for each selection string.
*/

class Parser {
  loos::Kernel krnl;
  ParserDriver driver;

public:
  explicit Parser(const string& s) : driver(s, krnl) { }    // This auto-parses the string...
  Parser() : driver(krnl) { }

  //! Parse the given command string after clearing the action stack
  void parse(const string& s) {krnl.clearActions();  driver.parse(s); }

  //! Return a ref to the compiled (hopefully) Kernel.
  loos::Kernel& kernel(void) { return(krnl); }
};





#endif

