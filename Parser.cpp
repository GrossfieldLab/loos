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
    if (t.type == Token::ID) {
      Action* r;
      if (t.datum == "name")
	r = new pushAtomName;
      else if (t.datum == "resname")
	r = new pushAtomResname;
      else if (t.datum == "segid")
	r = new pushAtomSegid;
      else {
	lex.push(t);
	return(0);
      }
	

      return(r);
    }

    lex.push(t);
    return(0);
  }

  Action* Parser::numId(void) {
    Token t = lex.pop();
    if (t.type == Token::ID) {
      Action* r;
      if (t.datum == "id")
	r = new pushAtomId;
      else if (t.datum == "resid")
	r = new pushAtomResid;
      else {
	lex.push(t);
	return(0);
      }
	

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
    return(r);
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


  Action* Parser::literal(void) {
    Action* r;

    r = number();
    if (!r)
      r = slit();

    return(r);
  }

  deque<Action*> Parser::stringexpr(void) {
    deque<Action*> acts;
    Action *a, *b, *c;

    a = slit();
    if (a) {
      b = stringop();
      if (b) {
	c = slit();
	if (c) {
	  acts.push_back(a);
	  acts.push_back(c);
	  acts.push_back(b);
	  return(acts);
	}
	delete b;
      }
      delete a;
    }

    acts.clear();
    return(acts);
  }


  deque<Action*> Parser::numexpr(void) {
    deque<Action*> acts;
    Action *a, *b, *c;

    a = number();
    if (a) {
      b = numop();
      if (b) {
	c = number();
	if (c) {
	  acts.push_back(a);
	  acts.push_back(c);
	  acts.push_back(b);
	  return(acts);
	}
	delete b;
      }
      delete a;
    }

    acts.clear();
    return(acts);
  }

  deque<Action*> Parser::subexpr(void) {
    deque<Action*> u;
    Action *r;

    u = numexpr();
    if (! u.empty())
      return(u);

    u = stringexpr();
    if (! u.empty())
      return(u);

    r = literal();
    if (r)
      u.push_back(r);

    return(u);
  }


  deque<Action*> Parser::expr(void) {
    deque<Action*> u, v, w;
    Action* r;
    Token t;

    r = unop();
    if (r) {
      u = expr();
      if (! u.empty()) {
	u.push_back(r);
	return(u);
      }
    }

    u = subexpr();
    if (! u.empty()) {
      r = binop();
      if (r) {
	v = expr();
	if (! v.empty()) {
	  w = u;
	  w.insert(w.end(), v.begin(), v.end());
	  w.push_back(r);
	  return(w);
	}
      } 
    }

    return(u);
  }


  void Parser::parse(const string s) {
    lex.tokenize(s);

    deque<Action*> acts = expr();
    if (acts.empty())
      throw(runtime_error("Parse failed"));

    deque<Action*>::iterator i;
    for (i=acts.begin(); i != acts.end(); i++)
      kernel->push(*i);
  }
  


};
