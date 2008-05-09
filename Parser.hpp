/*
  Parser.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Parser for atom selections...

  Current grammar is (roughly):


  <numeric> ::= [+-]\d+
  <alphabetic> ::= '[^']+'
  <id> ::= \w+
  <operator> ::= [&|!<>=]+

  <stringkey> ::= name | resname | segid
  <numkey> ::= id | resid

  <string> ::= <stringkey>|<alphabetic>
  <number> ::= <numkey>|<numeric>

  <stringop> ::= <|>|<=|>=|==
  <numberop> ::= <|>|<=|>=|==

  <stringexpr> ::= <string> <stringop> <string>
  <numexpr> ::= <number> <numberop> <number>

  <binop> ::= && | ||
  <unop> ::= !

  <literal> ::= <number> | <string>

  <subexpr> ::= <numexpr> | <stringexpr> | <literal>
  <expr> ::= <expr> <binop> <expr> | <unop> <expr> | <subexpr>

*/


#if !defined(PARSER_HPP)
#define PARSER_HPP

#include <iostream>
#include <string>
#include <vector>
#include <deque>


using namespace std;

#include "Tokenizer.hpp"
#include "Kernel.hpp"

namespace loos {

  class Parser {
    Kernel* kernel;
    Tokenizer lex;
    
  private:
    bool parsePrimNumeric(void);
    bool parsePrimId(void);
    bool parsePrimOperator(void);
    bool parsePrimString(void);

    bool parseStringKey(void);
    bool parseNumKey(void);

    bool parseString(void);
    bool parseNumber(void);

    bool parseStringExpr(void);
    bool parseNumExpr(void);

    bool parseBinop(void);
    bool parseUnop(void);

    bool parseLiteral(void);

    bool parseSubExpr(void);
    bool parseExpr(void);

  public:
    Parser(Kernel* k) : kernel(k) { }

    void parse(const string s);

  };

};
