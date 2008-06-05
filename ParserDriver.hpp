/*
  ParserDriver.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Driver for the Bison-generated parser...
*/


#if !defined(PARSERDRIVER_HPP)
#define PARSERDRIVER_HPP

#include <string>
#include <vector>

#include <FlexLexer.h>

#include "Kernel.hpp"

#include "grammar.hh"
#include "LoosLexer.hpp"

using namespace std;
using namespace loos;

//! Driver for the Bison parser (to encapsulate data)
//! Can parse from either stdin or a string.  Requires a Kernel for
//! storing the compiled Actions
struct ParserDriver {
  loos::parser *pparser;
  LoosLexer *lexer;
  Kernel& kern;
  istringstream *isp;
  

  //! For future parsing...
  ParserDriver(Kernel& k) : pparser(0), lexer(0), kern(k), isp(0) { }

  //! For parsing a string...
  ParserDriver(const string s, Kernel& k) : pparser(0), lexer(0), kern(k), isp(0) {
    if (isp)
      delete isp;
    isp = new istringstream(s);

    if (lexer)
      delete lexer;

    lexer = new LoosLexer(isp);
    parse();
  }

  ~ParserDriver() { delete pparser; delete lexer; delete isp; }

  //! Parse the passed string...
  /**
   *Note that it is up to the caller to reset the kernel if you don't
   *want to concatenate the commands...
   */
  void parse(const string& s) {
    if (isp)
      delete isp;
    isp = new istringstream(s);

    if (lexer)
      delete lexer;
    lexer = new LoosLexer(isp);

    parse();
  }

  //! Calls the Bison parser
  void parse(void) {
    if (!lexer)
      throw(runtime_error("Attempting to parse sans lexer"));

    if (!pparser)
      pparser = new loos::parser(*this);
    if (pparser->parse())
      throw(runtime_error("Parse error"));
  }

};





#endif
