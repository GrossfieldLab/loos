/* A Bison parser, made by GNU Bison 2.4.3.  */

/* Skeleton implementation for Bison LALR(1) parsers in C++
   
      Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 Free
   Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

// Take the name prefix into account.
#define yylex   looslex

/* First part of user declarations.  */

/* Line 311 of lalr1.cc  */
#line 6 "grammar.yy"


#include <ctype.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>

namespace loos {
	  struct ParserDriver;
};




/* Line 311 of lalr1.cc  */
#line 59 "grammar.cc"


#include "grammar.hh"

/* User implementation prologue.  */

/* Line 317 of lalr1.cc  */
#line 30 "grammar.yy"


#include "ParserDriver.hpp"
#include "LoosLexer.hpp"

#include "Kernel.hpp"

#undef yylex
#define yylex driver.lexer->looslex


namespace loos {
  void parse_error(const std::string&);
};



/* Line 317 of lalr1.cc  */
#line 86 "grammar.cc"

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* FIXME: INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#define YYUSE(e) ((void) (e))

/* Enable debugging if requested.  */
#if YYDEBUG

/* A pseudo ostream that takes yydebug_ into account.  */
# define YYCDEBUG if (yydebug_) (*yycdebug_)

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)	\
do {							\
  if (yydebug_)						\
    {							\
      *yycdebug_ << Title << ' ';			\
      yy_symbol_print_ ((Type), (Value), (Location));	\
      *yycdebug_ << std::endl;				\
    }							\
} while (false)

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug_)				\
    yy_reduce_print_ (Rule);		\
} while (false)

# define YY_STACK_PRINT()		\
do {					\
  if (yydebug_)				\
    yystack_print_ ();			\
} while (false)

#else /* !YYDEBUG */

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_REDUCE_PRINT(Rule)
# define YY_STACK_PRINT()

#endif /* !YYDEBUG */

#define yyerrok		(yyerrstatus_ = 0)
#define yyclearin	(yychar = yyempty_)

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)


namespace loos {

/* Line 380 of lalr1.cc  */
#line 152 "grammar.cc"
#if YYERROR_VERBOSE

  /* Return YYSTR after stripping away unnecessary quotes and
     backslashes, so that it's suitable for yyerror.  The heuristic is
     that double-quoting is unnecessary unless the string contains an
     apostrophe, a comma, or backslash (other than backslash-backslash).
     YYSTR is taken from yytname.  */
  std::string
  parser::yytnamerr_ (const char *yystr)
  {
    if (*yystr == '"')
      {
        std::string yyr = "";
        char const *yyp = yystr;

        for (;;)
          switch (*++yyp)
            {
            case '\'':
            case ',':
              goto do_not_strip_quotes;

            case '\\':
              if (*++yyp != '\\')
                goto do_not_strip_quotes;
              /* Fall through.  */
            default:
              yyr += *yyp;
              break;

            case '"':
              return yyr;
            }
      do_not_strip_quotes: ;
      }

    return yystr;
  }

#endif

  /// Build a parser object.
  parser::parser (ParserDriver& driver_yyarg)
    :
#if YYDEBUG
      yydebug_ (false),
      yycdebug_ (&std::cerr),
#endif
      driver (driver_yyarg)
  {
  }

  parser::~parser ()
  {
  }

#if YYDEBUG
  /*--------------------------------.
  | Print this symbol on YYOUTPUT.  |
  `--------------------------------*/

  inline void
  parser::yy_symbol_value_print_ (int yytype,
			   const semantic_type* yyvaluep, const location_type* yylocationp)
  {
    YYUSE (yylocationp);
    YYUSE (yyvaluep);
    switch (yytype)
      {
         default:
	  break;
      }
  }


  void
  parser::yy_symbol_print_ (int yytype,
			   const semantic_type* yyvaluep, const location_type* yylocationp)
  {
    *yycdebug_ << (yytype < yyntokens_ ? "token" : "nterm")
	       << ' ' << yytname_[yytype] << " ("
	       << *yylocationp << ": ";
    yy_symbol_value_print_ (yytype, yyvaluep, yylocationp);
    *yycdebug_ << ')';
  }
#endif

  void
  parser::yydestruct_ (const char* yymsg,
			   int yytype, semantic_type* yyvaluep, location_type* yylocationp)
  {
    YYUSE (yylocationp);
    YYUSE (yymsg);
    YYUSE (yyvaluep);

    YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

    switch (yytype)
      {
        case 4: /* "STRING" */

/* Line 480 of lalr1.cc  */
#line 74 "grammar.yy"
	{ delete (yyvaluep->sval); };

/* Line 480 of lalr1.cc  */
#line 259 "grammar.cc"
	break;
      case 5: /* "SKEY" */

/* Line 480 of lalr1.cc  */
#line 74 "grammar.yy"
	{ delete (yyvaluep->sval); };

/* Line 480 of lalr1.cc  */
#line 268 "grammar.cc"
	break;
      case 6: /* "NKEY" */

/* Line 480 of lalr1.cc  */
#line 74 "grammar.yy"
	{ delete (yyvaluep->sval); };

/* Line 480 of lalr1.cc  */
#line 277 "grammar.cc"
	break;
      case 30: /* "string" */

/* Line 480 of lalr1.cc  */
#line 74 "grammar.yy"
	{ delete (yyvaluep->sval); };

/* Line 480 of lalr1.cc  */
#line 286 "grammar.cc"
	break;
      case 31: /* "strval" */

/* Line 480 of lalr1.cc  */
#line 74 "grammar.yy"
	{ delete (yyvaluep->sval); };

/* Line 480 of lalr1.cc  */
#line 295 "grammar.cc"
	break;

	default:
	  break;
      }
  }

  void
  parser::yypop_ (unsigned int n)
  {
    yystate_stack_.pop (n);
    yysemantic_stack_.pop (n);
    yylocation_stack_.pop (n);
  }

#if YYDEBUG
  std::ostream&
  parser::debug_stream () const
  {
    return *yycdebug_;
  }

  void
  parser::set_debug_stream (std::ostream& o)
  {
    yycdebug_ = &o;
  }


  parser::debug_level_type
  parser::debug_level () const
  {
    return yydebug_;
  }

  void
  parser::set_debug_level (debug_level_type l)
  {
    yydebug_ = l;
  }
#endif

  int
  parser::parse ()
  {
    /// Lookahead and lookahead in internal form.
    int yychar = yyempty_;
    int yytoken = 0;

    /* State.  */
    int yyn;
    int yylen = 0;
    int yystate = 0;

    /* Error handling.  */
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// Semantic value of the lookahead.
    semantic_type yylval;
    /// Location of the lookahead.
    location_type yylloc;
    /// The locations where the error started and ended.
    location_type yyerror_range[3];

    /// $$.
    semantic_type yyval;
    /// @$.
    location_type yyloc;

    int yyresult;

    YYCDEBUG << "Starting parse" << std::endl;


    /* Initialize the stacks.  The initial state will be pushed in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystate_stack_ = state_stack_type (0);
    yysemantic_stack_ = semantic_stack_type (0);
    yylocation_stack_ = location_stack_type (0);
    yysemantic_stack_.push (yylval);
    yylocation_stack_.push (yylloc);

    /* New state.  */
  yynewstate:
    yystate_stack_.push (yystate);
    YYCDEBUG << "Entering state " << yystate << std::endl;

    /* Accept?  */
    if (yystate == yyfinal_)
      goto yyacceptlab;

    goto yybackup;

    /* Backup.  */
  yybackup:

    /* Try to take a decision without lookahead.  */
    yyn = yypact_[yystate];
    if (yyn == yypact_ninf_)
      goto yydefault;

    /* Read a lookahead token.  */
    if (yychar == yyempty_)
      {
	YYCDEBUG << "Reading a token: ";
	yychar = yylex (&yylval);
      }


    /* Convert token to internal form.  */
    if (yychar <= yyeof_)
      {
	yychar = yytoken = yyeof_;
	YYCDEBUG << "Now at end of input." << std::endl;
      }
    else
      {
	yytoken = yytranslate_ (yychar);
	YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
      }

    /* If the proper action on seeing token YYTOKEN is to reduce or to
       detect an error, take that action.  */
    yyn += yytoken;
    if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yytoken)
      goto yydefault;

    /* Reduce or error.  */
    yyn = yytable_[yyn];
    if (yyn <= 0)
      {
	if (yyn == 0 || yyn == yytable_ninf_)
	goto yyerrlab;
	yyn = -yyn;
	goto yyreduce;
      }

    /* Shift the lookahead token.  */
    YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

    /* Discard the token being shifted.  */
    yychar = yyempty_;

    yysemantic_stack_.push (yylval);
    yylocation_stack_.push (yylloc);

    /* Count tokens shifted since error; after three, turn off error
       status.  */
    if (yyerrstatus_)
      --yyerrstatus_;

    yystate = yyn;
    goto yynewstate;

  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
  `-----------------------------------------------------------*/
  yydefault:
    yyn = yydefact_[yystate];
    if (yyn == 0)
      goto yyerrlab;
    goto yyreduce;

  /*-----------------------------.
  | yyreduce -- Do a reduction.  |
  `-----------------------------*/
  yyreduce:
    yylen = yyr2_[yyn];
    /* If YYLEN is nonzero, implement the default value of the action:
       `$$ = $1'.  Otherwise, use the top of the stack.

       Otherwise, the following line sets YYVAL to garbage.
       This behavior is undocumented and Bison
       users should not rely upon it.  */
    if (yylen)
      yyval = yysemantic_stack_[yylen - 1];
    else
      yyval = yysemantic_stack_[0];

    {
      slice<location_type, location_stack_type> slice (yylocation_stack_, yylen);
      YYLLOC_DEFAULT (yyloc, slice, yylen);
    }
    YY_REDUCE_PRINT (yyn);
    switch (yyn)
      {
	  case 3:

/* Line 678 of lalr1.cc  */
#line 80 "grammar.yy"
    { driver.kern.push(new internal::logicalAnd); }
    break;

  case 4:

/* Line 678 of lalr1.cc  */
#line 81 "grammar.yy"
    { driver.kern.push(new internal::logicalOr); }
    break;

  case 6:

/* Line 678 of lalr1.cc  */
#line 86 "grammar.yy"
    { driver.kern.push(new internal::logicalNot); }
    break;

  case 7:

/* Line 678 of lalr1.cc  */
#line 87 "grammar.yy"
    { driver.kern.push(new internal::lessThan); }
    break;

  case 8:

/* Line 678 of lalr1.cc  */
#line 88 "grammar.yy"
    { driver.kern.push(new internal::lessThanEquals); }
    break;

  case 9:

/* Line 678 of lalr1.cc  */
#line 89 "grammar.yy"
    { driver.kern.push(new internal::greaterThanEquals); }
    break;

  case 10:

/* Line 678 of lalr1.cc  */
#line 90 "grammar.yy"
    { driver.kern.push(new internal::greaterThan); }
    break;

  case 11:

/* Line 678 of lalr1.cc  */
#line 91 "grammar.yy"
    { driver.kern.push(new internal::equals); }
    break;

  case 12:

/* Line 678 of lalr1.cc  */
#line 92 "grammar.yy"
    { driver.kern.push(new internal::equals); driver.kern.push(new internal::logicalNot); }
    break;

  case 13:

/* Line 678 of lalr1.cc  */
#line 93 "grammar.yy"
    { driver.kern.push(new internal::matchRegex(*((yysemantic_stack_[(3) - (3)].sval)))); }
    break;

  case 14:

/* Line 678 of lalr1.cc  */
#line 94 "grammar.yy"
    { driver.kern.push(new internal::logicalTrue); }
    break;

  case 15:

/* Line 678 of lalr1.cc  */
#line 95 "grammar.yy"
    { driver.kern.push(new internal::Hydrogen); }
    break;

  case 22:

/* Line 678 of lalr1.cc  */
#line 103 "grammar.yy"
    { driver.kern.push(new internal::extractNumber(*((yysemantic_stack_[(3) - (3)].sval)))); }
    break;

  case 23:

/* Line 678 of lalr1.cc  */
#line 105 "grammar.yy"
    { driver.kern.push(new internal::pushInt((yysemantic_stack_[(1) - (1)].ival))); }
    break;

  case 26:

/* Line 678 of lalr1.cc  */
#line 110 "grammar.yy"
    { (yyval.sval) = (yysemantic_stack_[(1) - (1)].sval); driver.kern.push(new internal::pushString(*((yysemantic_stack_[(1) - (1)].sval)))); }
    break;

  case 28:

/* Line 678 of lalr1.cc  */
#line 119 "grammar.yy"
    {
(yyval.sval) = (yysemantic_stack_[(1) - (1)].sval);
if (*((yysemantic_stack_[(1) - (1)].sval)) == "name")
   driver.kern.push(new internal::pushAtomName);
else if (*((yysemantic_stack_[(1) - (1)].sval)) == "resname")
   driver.kern.push(new internal::pushAtomResname);
else if (*((yysemantic_stack_[(1) - (1)].sval)) == "segid" || *((yysemantic_stack_[(1) - (1)].sval)) == "segname")
   driver.kern.push(new internal::pushAtomSegid);
else
   loos::parse_error("Unknown string keyword " + *((yysemantic_stack_[(1) - (1)].sval)));
}
    break;

  case 29:

/* Line 678 of lalr1.cc  */
#line 134 "grammar.yy"
    {
if (*((yysemantic_stack_[(1) - (1)].sval)) == "id")
   driver.kern.push(new internal::pushAtomId);
else if (*((yysemantic_stack_[(1) - (1)].sval)) == "resid")
   driver.kern.push(new internal::pushAtomResid);
else
   loos::parse_error("Unknown numeric keyword " + *((yysemantic_stack_[(1) - (1)].sval)));   
}
    break;



/* Line 678 of lalr1.cc  */
#line 624 "grammar.cc"
	default:
          break;
      }
    YY_SYMBOL_PRINT ("-> $$ =", yyr1_[yyn], &yyval, &yyloc);

    yypop_ (yylen);
    yylen = 0;
    YY_STACK_PRINT ();

    yysemantic_stack_.push (yyval);
    yylocation_stack_.push (yyloc);

    /* Shift the result of the reduction.  */
    yyn = yyr1_[yyn];
    yystate = yypgoto_[yyn - yyntokens_] + yystate_stack_[0];
    if (0 <= yystate && yystate <= yylast_
	&& yycheck_[yystate] == yystate_stack_[0])
      yystate = yytable_[yystate];
    else
      yystate = yydefgoto_[yyn - yyntokens_];
    goto yynewstate;

  /*------------------------------------.
  | yyerrlab -- here on detecting error |
  `------------------------------------*/
  yyerrlab:
    /* If not already recovering from an error, report this error.  */
    if (!yyerrstatus_)
      {
	++yynerrs_;
	error (yylloc, yysyntax_error_ (yystate));
      }

    yyerror_range[1] = yylloc;
    if (yyerrstatus_ == 3)
      {
	/* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

	if (yychar <= yyeof_)
	  {
	  /* Return failure if at end of input.  */
	  if (yychar == yyeof_)
	    YYABORT;
	  }
	else
	  {
	    yydestruct_ ("Error: discarding", yytoken, &yylval, &yylloc);
	    yychar = yyempty_;
	  }
      }

    /* Else will try to reuse lookahead token after shifting the error
       token.  */
    goto yyerrlab1;


  /*---------------------------------------------------.
  | yyerrorlab -- error raised explicitly by YYERROR.  |
  `---------------------------------------------------*/
  yyerrorlab:

    /* Pacify compilers like GCC when the user code never invokes
       YYERROR and the label yyerrorlab therefore never appears in user
       code.  */
    if (false)
      goto yyerrorlab;

    yyerror_range[1] = yylocation_stack_[yylen - 1];
    /* Do not reclaim the symbols of the rule which action triggered
       this YYERROR.  */
    yypop_ (yylen);
    yylen = 0;
    yystate = yystate_stack_[0];
    goto yyerrlab1;

  /*-------------------------------------------------------------.
  | yyerrlab1 -- common code for both syntax error and YYERROR.  |
  `-------------------------------------------------------------*/
  yyerrlab1:
    yyerrstatus_ = 3;	/* Each real token shifted decrements this.  */

    for (;;)
      {
	yyn = yypact_[yystate];
	if (yyn != yypact_ninf_)
	{
	  yyn += yyterror_;
	  if (0 <= yyn && yyn <= yylast_ && yycheck_[yyn] == yyterror_)
	    {
	      yyn = yytable_[yyn];
	      if (0 < yyn)
		break;
	    }
	}

	/* Pop the current state because it cannot handle the error token.  */
	if (yystate_stack_.height () == 1)
	YYABORT;

	yyerror_range[1] = yylocation_stack_[0];
	yydestruct_ ("Error: popping",
		     yystos_[yystate],
		     &yysemantic_stack_[0], &yylocation_stack_[0]);
	yypop_ ();
	yystate = yystate_stack_[0];
	YY_STACK_PRINT ();
      }

    yyerror_range[2] = yylloc;
    // Using YYLLOC is tempting, but would change the location of
    // the lookahead.  YYLOC is available though.
    YYLLOC_DEFAULT (yyloc, yyerror_range, 2);
    yysemantic_stack_.push (yylval);
    yylocation_stack_.push (yyloc);

    /* Shift the error token.  */
    YY_SYMBOL_PRINT ("Shifting", yystos_[yyn],
		     &yysemantic_stack_[0], &yylocation_stack_[0]);

    yystate = yyn;
    goto yynewstate;

    /* Accept.  */
  yyacceptlab:
    yyresult = 0;
    goto yyreturn;

    /* Abort.  */
  yyabortlab:
    yyresult = 1;
    goto yyreturn;

  yyreturn:
    if (yychar != yyempty_)
      yydestruct_ ("Cleanup: discarding lookahead", yytoken, &yylval, &yylloc);

    /* Do not reclaim the symbols of the rule which action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    while (yystate_stack_.height () != 1)
      {
	yydestruct_ ("Cleanup: popping",
		   yystos_[yystate_stack_[0]],
		   &yysemantic_stack_[0],
		   &yylocation_stack_[0]);
	yypop_ ();
      }

    return yyresult;
  }

  // Generate an error message.
  std::string
  parser::yysyntax_error_ (int yystate)
  {
    std::string res;
    YYUSE (yystate);
#if YYERROR_VERBOSE
    int yyn = yypact_[yystate];
    if (yypact_ninf_ < yyn && yyn <= yylast_)
      {
	/* Start YYX at -YYN if negative to avoid negative indexes in
	   YYCHECK.  */
	int yyxbegin = yyn < 0 ? -yyn : 0;

	/* Stay within bounds of both yycheck and yytname.  */
	int yychecklim = yylast_ - yyn + 1;
	int yyxend = yychecklim < yyntokens_ ? yychecklim : yyntokens_;
	int count = 0;
	for (int x = yyxbegin; x < yyxend; ++x)
	  if (yycheck_[x + yyn] == x && x != yyterror_)
	    ++count;

	// FIXME: This method of building the message is not compatible
	// with internationalization.  It should work like yacc.c does it.
	// That is, first build a string that looks like this:
	// "syntax error, unexpected %s or %s or %s"
	// Then, invoke YY_ on this string.
	// Finally, use the string as a format to output
	// yytname_[tok], etc.
	// Until this gets fixed, this message appears in English only.
	res = "syntax error, unexpected ";
	res += yytnamerr_ (yytname_[tok]);
	if (count < 5)
	  {
	    count = 0;
	    for (int x = yyxbegin; x < yyxend; ++x)
	      if (yycheck_[x + yyn] == x && x != yyterror_)
		{
		  res += (!count++) ? ", expecting " : " or ";
		  res += yytnamerr_ (yytname_[x]);
		}
	  }
      }
    else
#endif
      res = YY_("syntax error");
    return res;
  }


  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
  const signed char parser::yypact_ninf_ = -16;
  const signed char
  parser::yypact_[] =
  {
         7,   -16,   -16,   -16,   -16,   -16,   -16,     7,     7,    28,
     -16,    41,   -16,   -16,   -16,   -16,   -16,   -15,   -16,   -16,
      -5,    30,   -16,     7,     7,    27,    27,    27,    27,    27,
      27,     2,     2,   -16,   -16,   -16,   -16,    27,   -16,   -11,
     -16,   -16,   -16,   -16,   -16,   -16,   -16,   -16,    13
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned char
  parser::yydefact_[] =
  {
         0,    23,    26,    28,    29,    14,    15,     0,     0,     0,
       2,     0,    17,    19,    20,    18,    24,    25,    21,     6,
       0,     0,     1,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     5,    16,     3,     4,     0,     7,    25,
       8,     9,    10,    11,    12,    27,    13,    22,     0
  };

  /* YYPGOTO[NTERM-NUM].  */
  const signed char
  parser::yypgoto_[] =
  {
       -16,    40,    16,    -8,   -16,   -16,   -16,   -16,   -16,     3,
       1,   -16
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const signed char
  parser::yydefgoto_[] =
  {
        -1,     9,    10,    11,    12,    13,    14,    15,    16,    46,
      39,    18
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const unsigned char
  parser::yytable_[] =
  {
        21,    17,    31,    32,    23,    24,    45,    32,    17,    17,
       1,     2,     3,     4,     5,     6,    33,    38,    40,    41,
      42,    43,    44,    19,    17,    17,     7,     8,    22,    48,
       1,     2,     3,     4,    34,    47,     0,    23,    24,    35,
      36,    25,    26,    27,    28,    29,    30,    37,    20,     0,
       0,    34,    25,    26,    27,    28,    29,    30
  };

  /* YYCHECK.  */
  const signed char
  parser::yycheck_[] =
  {
         8,     0,    17,    18,     9,    10,     4,    18,     7,     8,
       3,     4,     5,     6,     7,     8,    21,    25,    26,    27,
      28,    29,    30,     7,    23,    24,    19,    20,     0,    37,
       3,     4,     5,     6,    21,    32,    -1,     9,    10,    23,
      24,    11,    12,    13,    14,    15,    16,    20,     8,    -1,
      -1,    21,    11,    12,    13,    14,    15,    16
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned char
  parser::yystos_[] =
  {
         0,     3,     4,     5,     6,     7,     8,    19,    20,    23,
      24,    25,    26,    27,    28,    29,    30,    32,    33,    24,
      23,    25,     0,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    21,    21,    24,    24,    20,    25,    32,
      25,    25,    25,    25,    25,     4,    31,    31,    25
  };

#if YYDEBUG
  /* TOKEN_NUMBER_[YYLEX-NUM] -- Internal symbol number corresponding
     to YYLEX-NUM.  */
  const unsigned short int
  parser::yytoken_number_[] =
  {
         0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
      40,    41
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned char
  parser::yyr1_[] =
  {
         0,    22,    23,    23,    23,    24,    24,    24,    24,    24,
      24,    24,    24,    24,    24,    24,    25,    25,    25,    25,
      26,    26,    27,    28,    29,    29,    30,    31,    32,    33
  };

  /* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
  const unsigned char
  parser::yyr2_[] =
  {
         0,     2,     1,     3,     3,     3,     2,     3,     3,     3,
       3,     3,     3,     3,     1,     1,     3,     1,     1,     1,
       1,     1,     3,     1,     1,     1,     1,     1,     1,     1
  };

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
  /* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
     First, the terminals, then, starting at \a yyntokens_, nonterminals.  */
  const char*
  const parser::yytname_[] =
  {
    "END", "error", "$undefined", "NUMBER", "STRING", "SKEY", "NKEY", "ALL",
  "HYDROGEN", "\"&&\"", "\"||\"", "\"<\"", "\"<=\"", "\">=\"", "\">\"",
  "\"==\"", "\"!=\"", "\"=~\"", "\"->\"", "\"!\"", "'('", "')'", "$accept",
  "expr", "rexpr", "value", "numeric", "numex", "number", "alpha",
  "string", "strval", "alphid", "numid", 0
  };
#endif

#if YYDEBUG
  /* YYRHS -- A `-1'-separated list of the rules' RHS.  */
  const parser::rhs_number_type
  parser::yyrhs_[] =
  {
        23,     0,    -1,    24,    -1,    23,     9,    24,    -1,    23,
      10,    24,    -1,    20,    23,    21,    -1,    19,    24,    -1,
      25,    11,    25,    -1,    25,    12,    25,    -1,    25,    13,
      25,    -1,    25,    14,    25,    -1,    25,    15,    25,    -1,
      25,    16,    25,    -1,    32,    17,    31,    -1,     7,    -1,
       8,    -1,    20,    25,    21,    -1,    26,    -1,    29,    -1,
      27,    -1,    28,    -1,    33,    -1,    32,    18,    31,    -1,
       3,    -1,    30,    -1,    32,    -1,     4,    -1,     4,    -1,
       5,    -1,     6,    -1
  };

  /* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
     YYRHS.  */
  const unsigned char
  parser::yyprhs_[] =
  {
         0,     0,     3,     5,     9,    13,    17,    20,    24,    28,
      32,    36,    40,    44,    48,    50,    52,    56,    58,    60,
      62,    64,    66,    70,    72,    74,    76,    78,    80,    82
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned char
  parser::yyrline_[] =
  {
         0,    79,    79,    80,    81,    85,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    99,    99,    99,    99,
     101,   101,   103,   105,   107,   107,   110,   112,   119,   134
  };

  // Print the state stack on the debug stream.
  void
  parser::yystack_print_ ()
  {
    *yycdebug_ << "Stack now";
    for (state_stack_type::const_iterator i = yystate_stack_.begin ();
	 i != yystate_stack_.end (); ++i)
      *yycdebug_ << ' ' << *i;
    *yycdebug_ << std::endl;
  }

  // Report on the debug stream that the rule \a yyrule is going to be reduced.
  void
  parser::yy_reduce_print_ (int yyrule)
  {
    unsigned int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    /* Print the symbols being reduced, and their result.  */
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1
	       << " (line " << yylno << "):" << std::endl;
    /* The symbols being reduced.  */
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT ("   $" << yyi + 1 << " =",
		       yyrhs_[yyprhs_[yyrule] + yyi],
		       &(yysemantic_stack_[(yynrhs) - (yyi + 1)]),
		       &(yylocation_stack_[(yynrhs) - (yyi + 1)]));
  }
#endif // YYDEBUG

  /* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
  parser::token_number_type
  parser::yytranslate_ (int t)
  {
    static
    const token_number_type
    translate_table[] =
    {
           0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      20,    21,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 57;
  const int parser::yynnts_ = 12;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 22;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 22;

  const unsigned int parser::yyuser_token_number_max_ = 274;
  const parser::token_number_type parser::yyundef_token_ = 2;


} // loos

/* Line 1054 of lalr1.cc  */
#line 1075 "grammar.cc"


/* Line 1056 of lalr1.cc  */
#line 145 "grammar.yy"



void loos::parser::error(const loos::location& loc, const std::string& s = "unknown error") {
  std::cerr << "***ERROR***  Bad selection syntax - " << s << std::endl;
}

void loos::parse_error(const std::string& s = "unknown error") {
  std::cerr << "***ERROR***  Bad selection syntax - " << s << std::endl;
}

