/*
  KernelValue.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

*/


#if !defined(TOKENIZER_HPP)
#define TOKENIZER_HPP



#include <iostream>
#include <string>
#include <stdexcept>
#include <deque>

#include <string.h>


using namespace std;


namespace loos {


  struct Token {
    enum TokenType { NONE, ID, NUMERIC, STRING, OPERATOR, LPAR, RPAR };

    TokenType type;
    string datum;

    Token() : type(NONE), datum("NONE") { }


    void setId(const string s) { datum = s; type = ID; }
    void setNumeric(const string s) { datum = s; type = NUMERIC; }
    void setString(const string s) { datum = s; type = STRING; }
    void setOperator(const string s) { datum = s; type = OPERATOR; }
    void setLpar(const string s) { datum = s; type = LPAR; }
    void setRpar(const string s) { datum = s; type = RPAR; }

    friend ostream& operator<<(ostream& os, const Token& t) {
      os << "<TOKEN TYPE='";
      switch(t.type) {
      case NONE: os << "NONE"; break;
      case ID: os << "ID"; break;
      case NUMERIC: os << "NUMERIC" ; break;
      case STRING: os << "STRING" ; break;
      case OPERATOR: os << "OPERATOR" ; break;
      case LPAR: os << "LPAR"; break;
      case RPAR: os << "RPAR"; break;
      default: throw(logic_error("Should never be here"));
      }

      os << "' DATA='" << t.datum << "' \\>";

      return(os);
    }


  };

  struct Tokens {
    Tokens() { }
    Tokens(const Tokens& t) : list(t.list) { }
    const Tokens& operator=(const Tokens& t) { list = t.list; return(*this); }
    
    deque<Token>& tokens(void) { return(list); }
    Token pop(void) {
      Token t = list.front();
      list.pop_front();
      return(t);
    }

    void push(const Token& t) { list.push_back(t); }

    friend ostream& operator<<(ostream& os, const Tokens& t) {
      deque<Token>::const_iterator i;
      for (i = t.list.begin(); i != t.list.end(); i++)
	os << *i << endl;
      return(os);
    }

    deque<Token> list;
  };

  Tokens tokenize(const string s);

};


#endif
