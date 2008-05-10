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
  <operator> ::= [&|!<>=]+

  <stringid> ::= name | resname | segid
  <numid> ::= id | resid

  <slit> ::= <stringkey>|<alphabetic>
  <number> ::= <numkey>|<numeric>

  <relop> ::= <|>|<=|>=
  <stringop> ::= <relop>|==  {regexp matching}
  <numberop> ::= <relop>|==

  <stringexpr> ::= <slit> <stringop> <slit>
  <numexpr> ::= <number> <numberop> <number>

  <binop> ::= && | ||
  <unop> ::= !

  <literal> ::= <number> | <slit>

  <subexpr> ::= <numexpr> | <stringexpr> | <literal>
  <expr> ::= <subexpr> <binop> <expr> | <unop> <expr> | (<expr>) | <subexpr>

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
    Action* numeric(void);
    Action* alphabetic(void);
    Action* stringId(void);
    Action* numId(void);
    Action* slit(void);
    Action* number(void);
    Action* relop(void);
    Action* stringop(void);
    Action* numop(void);
    Action* binop(void);
    Action* unop(void);
    Action* literal(void);

    deque<Action*> stringexpr(void);
    deque<Action*> numexpr(void);
    deque<Action*> subexpr(void);
    deque<Action*> expr(void);

  public:
    Parser(Kernel* k) : kernel(k) { }

    void parse(const string s);

  };

};


#endif
