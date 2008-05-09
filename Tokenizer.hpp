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
    enum TokenType { ID, NUMERIC, STRING, OPERATOR };

    TokenType type;
    string datum;


    void setId(const string s) { datum = s; type = ID; }
    void setNumeric(const string s) { datum = s; type = NUMERIC; }
    void setString(const string s) { datum = s; type = STRING; }
    void setOperator(const string s) { datum = s; type = OPERATOR; }

    friend ostream& operator<<(ostream& os, const Token& t) {
      os << "<TOKEN TYPE='";
      switch(t.type) {
      case ID: os << "ID"; break;
      case NUMERIC: os << "NUMERIC" ; break;
      case STRING: os << "STRING" ; break;
      case OPERATOR: os << "OPERATOR" ; break;
      default: throw(logic_error("Should never be here"));
      }

      os << "'>" << t.datum << "</TOKEN>";

      return(os);
    }


  };


  typedef deque<Token> Tokens;


  class Tokenizer {
    Tokens _tokens;
    string text;

  public:
    Tokenizer(const string s) : text(s) { tokenize(); }
    Tokenizer() { }

    void tokenize(void);
    void tokenize(const string s) { text = s; tokenize(); }

    Tokens& tokens(void) { return(_tokens); }


    // Note: we're popping/pushing to the front in this case...
    Token pop(void) {
      Token t = _tokens.front();
      _tokens.pop_front();
      return(t);
    }

    void push(const Token& t) {
      _tokens.push_front(t);
    }

    friend ostream& operator<<(ostream& os, const Tokenizer& t) {
      Tokens::const_iterator i;

      for (i=t._tokens.begin(); i != t._tokens.end(); i++) 
	os << *i << endl;

      return(os);
    }
  };


};



#endif
