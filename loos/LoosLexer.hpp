/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/





#if !defined(LOOS_LEXER_HPP)
#define LOOS_LEXER_HPP

#include <string>

// Putative patch for debian
#if !defined(EOF)
#define EOF    (-1)
#endif

// Flex's c++ support is daft...
#undef yyFlexLexer
#define yyFlexLexer LoosFlexLexer
#include <loos/FlexLexer.h>



#include "loos/grammar.hh"

// @cond TOOLS_INTERNAL

//! Extension of the Flex yyFlexLexer class.
class LoosLexer : public LoosFlexLexer {
public:
  LoosLexer() { }
  LoosLexer(std::istream* in) : LoosFlexLexer(in, 0) { }
  
  loos::parser::token_type looslex(loos::parser::semantic_type* yylval);
  
};

// @endcond

#endif
