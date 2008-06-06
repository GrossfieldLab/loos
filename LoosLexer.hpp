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





#if !defined(LOOSLEXER_HPP)
#define LOOSLEXER_HPP

#include <string>


// Flex's c++ support is daft...
#undef yyFlexLexer
#define yyFlexLexer LoosFlexLexer
#include <FlexLexer.h>



#include "grammar.hh"

using namespace std;



//! Extension of the Flex yyFlexLexer class.
class LoosLexer : public LoosFlexLexer {
public:
  LoosLexer() { }
  LoosLexer(istream* in) : LoosFlexLexer(in, 0) { }

  loos::parser::token_type looslex(loos::parser::semantic_type* yylval);

};


#endif
