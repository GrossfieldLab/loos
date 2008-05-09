/*
  Parser.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Parser for atom selections...
*/

#include <sstream>

#include "Parser.hpp"

using namespace std;


namespace loos {

  Action* Parser::numeric(void) {
    Token t = lex.pop();
    if (t.type == Token::NUMERIC) {
      int i;
      if (! (stringstream(t.datum) >> i)) 
	throw(runtime_error("Cannot convert " + t.datum + " to a numeric"));
      Action *r = new pushInt(i);
      return(r);
    }

    lex.push(t);
    return(0);
  }

  Action* Parser::alphabetic(void) {
    Token t = lex.pop();
    if (t.type == Token::STRING) {
      Action* r = new pushString(t.datum);
      return(r);
    }

    lex.push(t);
    return(0);
  }

  Action* Parser::stringId(void) {
    Token t = lex.pop();
    if (t.type == Token::STRING) {
      Action* r;
      if (t.datum == "name")
	r = new pushAtomName;
      else if (t.datum == "resname")
	r = new pushAtomResname;
      else if (t.datum == "segid")
	r = new pushAtomSegid;
      else
	throw(runtime_error("Unknown identifier " + t.datum));

      return(r);
    }

    lex.push(t);
    return(0);
  }

  Action* Parser::numId(void) {
    Token t = lex.pop();
    if (t.type == Token::STRING) {
      Action* r;
      if (t.datum == "id")
	r = new pushAtomId;
      else if (t.datum == "resid")
	r = new pushAtomResid;
      else
	throw(runtime_error("Unknown identifier " + t.datum));

      return(r);
    }

    lex.push(t);
    return(0);
  }


  Action* Parser::slit(void) {
    Action* r;

    r = stringId();
    if (!r)
      r = alphabetic();
    return(r);
  }

  Action* Parser::number(void) {
    Action* r;

    r = numId();
    if (!r)
      r = numeric();
    return(r);
  }

  Action* Parser::relop(void) {
    Action* r;
    Token t = lex.pop();
    if (t.type == Token::OPERATOR) {
      if (t.datum == "<")
	r = new lessThan;
      else if (t.datum == "<=")
	r = new lessThanEquals;
      else if (t.datum == ">=")
	r = new greaterThanEquals;
      else if (t.datum == ">")
	r = new greaterThan;
      else {
	lex.push(t);
	return(0);
      }

      return(r);
    }

    lex.push(t);
    return(r);
  }

  
  Action* Parser::stringop(void) {
    Action* r;

    r = relop();
    if (!r) {
      Token t = lex.pop();
      if (t.type == Token::OPERATOR) {
	if (t.datum == "==") {
	  r = new matchStringAsRegex;
	  return(r);
	}
      }
      lex.push(t);
    }

    return(0);
  }

  Action* Parser::numop(void) {
    Action* r;
    
    r = relop();
    if (!r) {
      Token t = lex.pop();
      if (t.type == Token::OPERATOR) {
	if (t.datum == "==") {
	  r = new equals;
	  return(r);
	}
      }
      lex.push(t);
    }
    return(0);
  }

  Action* Parser::binop(void) {
    Token t = lex.pop();

    if (t.type == Token::OPERATOR) {
      if (t.datum == "&&")
	return(new logicalAnd);
      else if (t.datum == "||")
	return(new logicalOr);
    }

    lex.push(t);
    return(0);
  }

  Action* Parser::unop(void) {
    Token t = lex.pop();
    
    if (t.type == Token::OPERATOR)
      if (t.datum == "!")
	return(new logicalNot);
    lex.push(t);
    return(0);
  }


  void Parser::parse(const string s) {
    lex.tokenize(s);

  }
  


};
