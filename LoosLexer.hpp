#if !defined(MYLEX_HPP)
#define MYLEX_HPP

#include <string>



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
