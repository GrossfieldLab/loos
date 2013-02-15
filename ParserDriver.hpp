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




#if !defined(LOOS_PARSERDRIVER_HPP)
#define LOOS_PARSERDRIVER_HPP

#include <string>
#include <vector>

#include <FlexLexer.h>

#include "exceptions.hpp"
#include "Kernel.hpp"

namespace loos { struct ParserDriver; }

#include "grammar.hh"
#include "LoosLexer.hpp"


namespace loos {

  //! Driver for the Bison parser (to encapsulate data)
  //! Can parse from either stdin or a string.  Requires a Kernel for
  //! storing the compiled Actions
  struct ParserDriver {
    parser *pparser;
    LoosLexer *lexer;
    Kernel& kern;
    std::istringstream *isp;
  

    //! For future parsing...
    explicit ParserDriver(Kernel& k) : pparser(0), lexer(0), kern(k), isp(0) { }

    //! For parsing a string...
    ParserDriver(const std::string s, Kernel& k) : pparser(0), lexer(0), kern(k), isp(0) {
      if (isp)
        delete isp;
      isp = new std::istringstream(s);

      if (lexer)
        delete lexer;

      lexer = new LoosLexer(isp);
      parse();
    }

    ~ParserDriver() { delete pparser; delete lexer; delete isp; }

    //! Parse the passed string...
    /**
     *Note that it is up to the caller to reset the kernel if you don't
     *want to concatenate the commands...
     */
    void parse(const std::string& s) {
      if (isp)
        delete isp;
      isp = new std::istringstream(s);

      if (lexer)
        delete lexer;
      lexer = new LoosLexer(isp);

      parse();
    }

    //! Calls the Bison parser
    void parse(void) {
      if (!lexer)
        throw(std::runtime_error("Attempting to parse sans lexer"));

      if (!pparser)
        pparser = new parser(*this);
      if (pparser->parse())
        throw(ParseError("Parse error"));
    }

  };


}


#endif
