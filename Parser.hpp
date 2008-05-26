/*
  Parser.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Parser for atom selections using Bison/Flex

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
 *  name, resname, and segid, but NOT an id nor a resid.
 *
 *  String equality matches the entire string.  If you want to match a
 *  subset, you should use the '=~' operator.  In other words,
 *  \verbatim
 *    "CA" == "C"   -> false
 *    "C"  == "C"   -> true
 *    "CA" =~ "C"   -> true
 *  \endverbatim
 *
 *  Finally, the standard precedence and associativity that apply in
 *  C++ apply here.  Expressions are evaluated left to right and
 *  parenthesis may be used to alter precedence/evaluation order.
 *  Unlike C/C++ however, the logical operators do not short-circuit.
 *
 *  If there is a syntax error in the selection string, then a
 *  runtime_error() is thrown.
*/

class Parser {
  loos::Kernel krnl;
  ParserDriver driver;

public:
  Parser(const string& s) : driver(s, krnl) { }

  //! Return a ref to the compiled (hopefully) Kernel.
  Kernel& kernel(void) { return(krnl); }
};





#endif

