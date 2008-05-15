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
#include <string>
#include <vector>
#include <deque>

#include <boost/tuple/tuple.hpp>

using namespace std;

#include "Tokenizer.hpp"
#include "Kernel.hpp"

namespace loos {

  static const char* ptok_type_bindings[18] = {"none", "num", "alpha", "name", "id", "resid", "resname", "segid", "lt", "lte", "gte", "gt", "eq", "neq", "lnot", "land", "lor", "regex"};

  struct ptok {
    enum ptok_type { none = 0, num, alpha, name, id, resid, resname, segid, lt, lte, gte, gt, eq, neq, lnot, land, lor, regex };

    ptok() : type(none), val("") { }
    ptok(const ptok& p) : type(p.type), val(p.val) { }
    const ptok& operator=(const ptok& p) { type=p.type; val=p.val; return(*this); }

    friend ostream& operator<<(ostream& os, const ptok& p) {
      os << "ptok(" << p.type << ":" << ptok_type_bindings[p.type] << ") = '" << p.val << "'";
      return(os);
    }

    ptok_type type;
    string val;
  };

  typedef vector<ptok> ptoks;
  typedef boost::tuple<ptoks, Tokens> pares;

  namespace parser {
    void prptoks(const ptoks&);

    pares alpha(const pares);
    pares number(const pares);
    pares numid(const pares);
    pares alphid(const pares);

    pares numeric(const pares);
    pares alphabetic(const pares);
    pares literal(const pares);

    pares unop(const pares);
    pares binop(const pares);
    pares relop(const pares);
    pares regexp(const pares);

    pares rexp(const pares);
    pares lexp(const pares);
    pares sexp(const pares);
    pares expr(const pares);

    bool translate(const ptoks, Kernel&);
    bool parse(const string, Kernel&);
  };


};
  
  



#endif

