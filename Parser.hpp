/*
  Parser.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Parser for atom selections...

  Current grammar is (roughly):


  expr       ::= sexp | '(' sexp ')' | unop sexp | sexp binop sexp
  sexp       ::= '(' lexp ')' | lexp
  lexp       ::= rexp | unop rexp | rexp binop rexp
  rexp       ::= regexp | lit relop lit
  relop      ::= > >= <= == !=
  regexp     ::= alphid '=~' alpha
  binop      ::= && ||
  unop       ::= !
  lit        ::= alphabetic | numeric
  alphabetic ::= alpha | alphid
  alphid     ::= name | resname | segid
  numeric    ::= number | numid
  numid      ::= id | resid

*/


#if !defined(PARSER_HPP)
#define PARSER_HPP

#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

#include <Kernel.hpp>
#include <ParserDriver.hpp>



class Parser {
  loos::Kernel krnel;
  ParserDriver driver;

public:
  Parser(const string& s) : driver(s, krnel) { driver.parse(); }

  bool operator()(const pAtom& pa) {
    krnel.execute(pa);
    if (krnel.stack().size() != 1)
      throw(runtime_error("Execution error - unexpected values on stack"));

    loos::Value results = krnel.stack().pop();
    if (results.type != loos::Value::INT)
      throw(runtime_error("Execution error - unexpected value on top of stack"));

    return(results.itg);
  }

  Kernel& kernel(void) { return(krnel); }
};



#endif

