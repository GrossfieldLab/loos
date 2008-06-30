%{
#include <string>

#define yyterminate()     return token::END

#include "grammar.hh"
#include "LoosLexer.hpp"

typedef loos::parser::token token;
typedef loos::parser::token_type token_type;

#define YY_NO_UNISTD_H

%}


%option noyywrap nounput c++

ws       [ \t]+
number   [0-9]+
ID       [a-zA-Z][a-zA-Z0-9]+

%%
[\n+]    /* Swallow newlines? */

#.+
{number}                { yylval->ival = atoi(yytext); return(token::NUMBER); }
"<"                     { return(token::LT); }
"<="                    { return(token::LTE); }
">="                    { return(token::GTE); }
">"                     { return(token::GT); }
"=="                    { return(token::EQ); }
"!="                    { return(token::NE); }
"=~"                    { return(token::REGEX); }

"&&"                    { return(token::AND); }
"||"                    { return(token::OR); }

"->"                 { return(token::NEKEY); }

all                  { return(token::ALL); }
name|resname|segid   { yylval->sval = new string(yytext, yyleng); return(token::SKEY); }
id|resid             { yylval->sval = new string(yytext, yyleng); return(token::NKEY); }

\"|\'                {                /* Special handling for strings... */
 string delim(yytext, yyleng);
 int c;
 string text;

 while ((c = yyinput()) > 0) {
  if (c == '\n' || c == delim[0])
     break;
   text += c;
 }

 if (c < 0)
   return(token::END);
 yylval->sval = new string(text);
 return(token::STRING);
}

                       

{ws}                 /* Swallow up remaining white-space */

.                    { return(static_cast<token_type>(*yytext)); }

%%

int yyFlexLexer::yylex() {
    std::cerr << "Should never be here!\n";
    exit(-1);
}


/* OS X Flex install requires this (2.5.35) to co-exist with Linux 2.5.33 install */

#if !defined(yywrap)
int LoosFlexLexer::yywrap() { return(1); }
#endif
