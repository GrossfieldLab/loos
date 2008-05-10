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

  Parser::Parsed Parser::number(Parsed pre) {
    Actions acts = boost::get<0>(pre);
    Tokenizer toks = boost::get<1>(pre);

    Token t = toks.pop();
    if (t.type == Token::NUMERIC) {
      int i;
      if (! (stringstream(t.datum) >> i)) 
	throw(runtime_error("Cannot convert " + t.datum + " to a numeric"));
      acts.push_back(new pushInt(i));
    } else
      toks.restore();

    Parsed post(acts, toks);
    return(post);
  }

  Parser::Parsed Parser::numid(Parsed pre) {
    Actions acts = boost::get<0>(pre);
    Tokenizer toks = boost::get<1>(pre);

    Token t = toks.pop();
    if (t.type == Token::ID) {
      Action* r;
      if (t.datum == "id")
	r = new pushAtomId;
      else if (t.datum == "resid")
	r = new pushAtomResid;
      else
	return(pre);
      acts.push_back(r);
      } else
	toks.restore();

    Parsed post(acts, toks);
    return(post);
  }
    
  Parser::Parsed Parser::alpha(Parsed pre) {
    Actions acts = boost::get<0>(pre);
    Tokenizer toks = boost::get<1>(pre);

    Token t = toks.pop();
    if (t.type == Token::STRING) {
      acts.push_back(new pushString(t.datum));
    } else
      toks.restore();

    Parsed post(acts, toks);
    return(post);
  }

  Parser::Parsed Parser::alphid(Parsed pre) {
    Actions acts = boost::get<0>(pre);
    Tokenizer toks = boost::get<1>(pre);

    Token t = toks.pop();
    if (t.type == Token::ID) {
      Action* r;
      if (t.datum == "name")
	r = new pushAtomName;
      else if (t.datum == "resname")
	r = new pushAtomResname;
      else if (t.datum == "segid")
	r = new pushAtomSegid;
      else
	return(pre);
      acts.push_back(r);
      } else
	toks.restore();

    Parsed post(acts, toks);
    return(post);
  }


  int Parser::depth(const Parsed& p) {
    Actions acts = boost::get<0>(p);
    return(acts.size());
  }
    


  Parser::Parsed Parser::numeric(Parsed pre) {
    Parsed opt1 = number(pre);
    Parsed opt2 = numid(pre);

    return( (depth(opt1) >= depth(opt2)) ? opt1 : opt2 );
  }




  Parser::Parsed Parser::alphabetic(Parsed pre) {
    Parsed opt1 = alpha(pre);
    Parsed opt2 = alphid(pre);

    return( (depth(opt1) >= depth(opt2)) ? opt1 : opt2 );
  }

  Parser::Parsed Parser::lit(Parsed pre) {
    Parsed opt1 = numeric(pre);
    Parsed opt2 = alphabetic(pre);

    return( (depth(opt1) >= depth(opt2)) ? opt1 : opt2 );
  }

  Parser::Parsed Parser::regexpr(Parsed pre) {
    Parsed op1 = alphid(pre);

    Actions acts = boost::get<0>(op1);
    Tokenizer toks = boost::get<1>(op1);

    Token t = toks.pop();
    if (t.type == Token::OPERATOR && t.datum == "=~") {
      t = toks.pop();
      if (t.type == Token::STRING) {
	acts.push_back(new matchRegex(t.datum));
	Parsed post(acts, toks);
	return(post);
      }
    }

    return(pre);
  }
  



  void Parser::parse(const string s) {
    cout << "Hello world\n";

  }

};


