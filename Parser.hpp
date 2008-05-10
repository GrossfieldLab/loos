/*
  Parser.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Parser for atom selections...

  Current grammar is (roughly):


       <exp> ::= <sexp> | <unop> <sexp> | <sexp> <binop> <sexp>
      <sexp> ::= <lexp> | '(' <lexp> ')'
      <lexp> ::= <rexp> | <unop> <rexp> | <rexp> <binop> <rexp>
      <rexp> ::= <lit> <relop> <lit> | <regexpr>
   <regexpr> ::= <alphid> '=~' <alpha>

     <relop> ::= '<' | '>' | '<=' | '>=' | '=='
     <binop> ::= '&&' | '||'
      <unop> ::= '!'
       <lit> ::= <numeric>|<alphabetic>
   <numeric> ::= <number>|<numid>
<alphabetic> ::= <alpha>|<alphid>

*/


#if !defined(PARSER_HPP)
#define PARSER_HPP

#include <iostream>
#include <string>
#include <vector>
#include <deque>

#include <boost/tuple/tuple.hpp>

using namespace std;

#include "Tokenizer.hpp"
#include "Kernel.hpp"

namespace loos {

  class Parser {

    typedef deque<Action*> Actions;
    typedef boost::tuple<Actions, Tokenizer> Parsed;

    Kernel* kernel;
    
  private:
    int depth(const Parsed& p);

    Parsed alpha(Parsed);
    Parsed alphid(Parsed);
    Parsed number(Parsed);
    Parsed numid(Parsed);
    Parsed numeric(Parsed);
    Parsed alphabetic(Parsed);
    Parsed lit(Parsed);
    Parsed regexpr(Parsed);


  public:
    Parser(Kernel* k) : kernel(k) { }

    void parse(const string s);

  };

};


#endif
