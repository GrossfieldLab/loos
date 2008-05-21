%{
#include <string>

#define yyterminate()     return token::END

#include "grammar.tab.hh"
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

name|resname|segid   { yylval->sval = new string(yytext, yyleng); return(token::SKEY); }
id|resid             { yylval->sval = new string(yytext, yyleng); return(token::NKEY); }
\"|\'                {
 string delim(yytext, yyleng);
 int c;
 string text;

 while ((c = yyinput())) {
   if (c == '\\')
     c = yyinput();
   else
     if (c == '\n' || c == delim[0])
       break;
   text += c;
 }
 yylval->sval = new string(text);
 return(token::STRING);
}

                       

{ws}

.                    { return(static_cast<token_type>(*yytext)); }

%%

int yyFlexLexer::yylex() {
    std::cerr << "Should never be here!\n";
    exit(-1);
}
