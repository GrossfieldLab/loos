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


//using namespace std;
//using namespace loos;

class ParserDriver;


#define YY_DECL loos::parser::token_type LoosLexer::looslex(loos::parser::semantic_type* yylval)


%}

%parse-param { ParserDriver& driver }

%union
{
	std::string *sval;
	int ival;
};

%{

#include "ParserDriver.hpp"
#include "LoosLexer.hpp"

#include "Kernel.hpp"

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

%left "&&" "||"
%left "<" "<=" ">=" ">" "==" "!=" "=~"

%type <sval> string alphid strval

%destructor { delete $$; } STRING string strval SKEY NKEY

%%


expr  : rexpr
      | expr "&&" rexpr   { driver.kern.push(new logicalAnd); }
      | expr "||" rexpr   { driver.kern.push(new logicalOr); }
      ;


rexpr : '(' expr ')'
      | '!' rexpr              { driver.kern.push(new logicalNot); }
      | value "<" value { driver.kern.push(new lessThan); }
      | value "<=" value { driver.kern.push(new lessThanEquals); }
      | value ">=" value { driver.kern.push(new greaterThanEquals); }
      | value ">" value { driver.kern.push(new greaterThan); }
      | value "==" value { driver.kern.push(new equals); }
      | value "!=" value { driver.kern.push(new equals); driver.kern.push(new logicalNot); }
      | alphid "=~" strval { driver.kern.push(new matchRegex(*($3))); }
      | ALL { driver.kern.push(new logicalTrue); }
      | HYDROGEN { driver.kern.push(new Hydrogen); }
      ;


value : '(' value ')' | numeric | alpha | numex;

numeric : number | numid;

numex : alphid "->" strval { driver.kern.push(new extractNumber(*($3))); } ;

number  : NUMBER        { driver.kern.push(new pushInt($1)); }  ;

alpha   : string | alphid
        ;

string  : STRING        { $$ = $1; driver.kern.push(new pushString(*($1))); } ;

strval  : STRING ;      /* Non-pushing string so we can compiled it */
                        /* into the regular expression operator */


/* Since we only have a few keywords, we just manually extract them
rather than use a table... */

alphid  : SKEY          {
$$ = $1;
if (*($1) == "name")
   driver.kern.push(new pushAtomName);
else if (*($1) == "resname")
   driver.kern.push(new pushAtomResname);
else if (*($1) == "segid" || *($1) == "segname")
   driver.kern.push(new pushAtomSegid);
else
   loos::parse_error("Unknown string keyword " + *($1));
}
        ;      



numid   : NKEY          {
if (*($1) == "id")
   driver.kern.push(new pushAtomId);
else if (*($1) == "resid")
   driver.kern.push(new pushAtomResid);
else
   loos::parse_error("Unknown numeric keyword " + *($1));   
}
        ;


%%


void loos::parser::error(const loos::location& loc, const std::string& s = "unknown error") {
  cerr << "***ERROR***  Bad selection syntax - " << s << endl;
}

void loos::parse_error(const std::string& s = "unknown error") {
  cerr << "***ERROR***  Bad selection syntax - " << s << endl;
}
