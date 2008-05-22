/*
  LoosLexer.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Extension of yyFlexLexer for parsing atom selections...
*/



#if !defined(LOOSLEXER_HPP)
#define LOOSLEXER_HPP

#include <string>


// Flex's c++ support is daft...
#undef yyFlexLexer
#define yyFlexLexer LoosFlexLexer
#include <FlexLexer.h>



#include "grammar.hh"

using namespace std;




class LoosLexer : public LoosFlexLexer {
public:
  LoosLexer() { }
  LoosLexer(istream* in) : LoosFlexLexer(in, 0) { }

  loos::parser::token_type looslex(loos::parser::semantic_type* yylval);

};


#endif
