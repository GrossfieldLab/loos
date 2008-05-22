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


struct ParserDriver {
  loos::parser *pparser;
  LoosLexer *lexer;
  Kernel& kern;
  istringstream *isp;
  
  vector<string> cmds;


  // For parsing from stdin
  ParserDriver(Kernel& k) : pparser(0), lexer(new LoosLexer), kern(k), isp(0) { parse(); }

  // For parsing a string...
  ParserDriver(const string s, Kernel& k) : pparser(0), lexer(0), kern(k), isp(0) {
    isp = new istringstream(s);
    lexer = new LoosLexer(isp);
    parse();
  }

  ~ParserDriver() { delete pparser; delete lexer; delete isp; }

  
  void parse(void) {
    if (!pparser)
      pparser = new loos::parser(*this);
    if (pparser->parse())
      throw(runtime_error("Parse error"));
  }

};





#endif
