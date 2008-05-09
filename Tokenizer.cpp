/*
  KernelValue.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

*/

#include <ctype.h>



#include "Tokenizer.hpp"


namespace loos {

  void Tokenizer::tokenize(void) {
    int i = 0;
    int state = 0;
    char c;
    string t;
    Token token;

    while (i < text.size()) {
      c = text[i];

      if (state == 0) {

	if (isalpha(c)) {
	  state = 1;
	  t = c;
	} else if (c == '-' || c == '+' || isdigit(c)) {
	  state = 2;
	  t = c;
	} else if (c == '\'') {
	  state = 3;
	  t = string("");
	} else if (c == '&' || c == '|' || c == '!' || c == '=' || c == '<' || c == '>') {
	  state = 4;
	  t = c;
	}

      } else if (state == 1) {

	if (!isalpha(c)) {
	  state = 0;
	  token.setId(t);
	  _tokens.push_back(token);
	  --i;
	} else 
	  t += c;

      } else if (state == 2) {

	if (!isdigit(c)) {
	  state = 0;
	  token.setNumeric(t);
	  _tokens.push_back(token);
	  --i;
	} else 
	  t += c;

      } else if (state == 3) {

	if (c == '\\') {
	  c = text[++i];
	  t += c;
	} else if (c == '\'') {
	  state = 0;
	  token.setString(t);
	  _tokens.push_back(token);
	} else
	  t += c;

      } else if (state == 4) {

	if (!(c == '&' || c == '|' || c == '!' || c == '=' || c == '<' || c == '>')) {
	  state = 0;
	  
	  // Validate operator...
	  if (!(t == "==" ||
		t == "&&" ||
		t == "||" ||
		t == "!" ||
		t == "<" ||
		t == ">" ||
		t == ">=" ||
		t == "<="))
	    throw(runtime_error("Unidentifed operator: " + t));

	  token.setOperator(t);
	  _tokens.push_back(token);
	  --i;
	} else
	  t += c;

      }


      ++i;
    }

    switch(state) {
    case 0: break;
    case 1: token.setId(t); _tokens.push_back(token); break;
    case 2: token.setNumeric(t); _tokens.push_back(token); break;
    case 3: token.setString(t); _tokens.push_back(token); break;
    case 4: token.setOperator(t); _tokens.push_back(token); break;
    default:
      throw(logic_error("Invalid state in tokenization"));
    }

  }

};
