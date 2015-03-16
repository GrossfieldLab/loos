%skeleton "lalr1.cc"
%defines
%name-prefix="loos"


%{

#include <ctype.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>

namespace loos {
	  struct ParserDriver;
};


%}

%parse-param { ParserDriver& driver }

%union
{
	std::string *sval;
	int ival;
};

%{

#include "loos/ParserDriver.hpp"
#include "loos/LoosLexer.hpp"

#include "loos/Kernel.hpp"

#undef yylex
#define yylex driver.lexer->looslex


namespace loos {
  void parse_error(const std::string&);
};

%}


%token		END	0
%token <ival>	NUMBER
%token <sval>	STRING
%token <sval>   SKEY
%token <sval>   NKEY
%token <sval>   ALL
%token <sval>   HYDROGEN
%token <sval>   BACKBONE


%token AND "&&"
%token OR  "||"
%token LT  "<"
%token LTE "<="
%token GTE ">="
%token GT  ">"
%token EQ  "=="
%token NE  "!="
%token REGEX "=~"
%token NEKEY "->"
%token NOT "!"

%left "&&" "||"
%left "<" "<=" ">=" ">" "==" "!=" "=~"

%type <sval> string alphid strval

%destructor { delete $$; } STRING string strval SKEY NKEY

%%


expr  : rexpr
      | expr "&&" rexpr   { driver.kern.push(new internal::logicalAnd); }
      | expr "||" rexpr   { driver.kern.push(new internal::logicalOr); }
      ;


rexpr : '(' expr ')'
      | NOT rexpr              { driver.kern.push(new internal::logicalNot); }
      | value "<" value { driver.kern.push(new internal::lessThan); }
      | value "<=" value { driver.kern.push(new internal::lessThanEquals); }
      | value ">=" value { driver.kern.push(new internal::greaterThanEquals); }
      | value ">" value { driver.kern.push(new internal::greaterThan); }
      | value "==" value { driver.kern.push(new internal::equals); }
      | value "!=" value { driver.kern.push(new internal::equals); driver.kern.push(new internal::logicalNot); }
      | alphid "=~" strval { driver.kern.push(new internal::matchRegex(*($3))); }
      | ALL { driver.kern.push(new internal::logicalTrue); }
      | HYDROGEN { driver.kern.push(new internal::Hydrogen); }
      | BACKBONE { driver.kern.push(new internal::Backbone); }
      ;


value : '(' value ')' | numeric | alpha | numex;

numeric : number | numid;

numex : alphid "->" strval { driver.kern.push(new internal::extractNumber(*($3))); } ;

number  : NUMBER        { driver.kern.push(new internal::pushInt($1)); }  ;

alpha   : string | alphid
        ;

string  : STRING        { $$ = $1; driver.kern.push(new internal::pushString(*($1))); } ;

strval  : STRING ;      /* Non-pushing string so we can compiled it */
                        /* into the regular expression operator */


/* Since we only have a few keywords, we just manually extract them
rather than use a table... */

alphid  : SKEY          {
$$ = $1;
if (*($1) == "name")
   driver.kern.push(new internal::pushAtomName);
else if (*($1) == "resname")
   driver.kern.push(new internal::pushAtomResname);
else if (*($1) == "segid" || *($1) == "segname")
   driver.kern.push(new internal::pushAtomSegid);
else if (*($1) == "chainid")
   driver.kern.push(new internal::pushAtomChainId);
else
   loos::parse_error("Unknown string keyword " + *($1));
}
        ;      



numid   : NKEY          {
if (*($1) == "id")
   driver.kern.push(new internal::pushAtomId);
else if (*($1) == "resid")
   driver.kern.push(new internal::pushAtomResid);
else
   loos::parse_error("Unknown numeric keyword " + *($1));   
}
        ;


%%


void loos::parser::error(const loos::location& loc, const std::string& s = "unknown error") {
  std::cerr << "***ERROR***  Bad selection syntax - " << s << std::endl;
}

void loos::parse_error(const std::string& s = "unknown error") {
  std::cerr << "***ERROR***  Bad selection syntax - " << s << std::endl;
}
