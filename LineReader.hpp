/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2012, Tod D. Romo, Alan Grossfield
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


#if !defined(LOOS_LINE_READER_HPP)
#define LOOS_LINE_READER_HPP


#include <iostream>
#include <string>
#include <stdexcept>
#include <list>
#include <loos_defs.hpp>


namespace loos {


  class LineReader {
  public:
    LineReader(std::istream& is) : _is(is), _lineno(0) { }
    LineReader(std::istream& is, std::string& name) : _is(is), _name(name) { }

    virtual bool getNext();

    virtual void push_back(const std::string& s);

    virtual bool eof() const;
    virtual bool fail() const;
    virtual bool good() const;
    
    virtual std::string line() const;
    
    virtual uint lineNumber() const;

  protected:

    virtual void checkState() const;

    virtual void stripComment(std::string& s) const;

    virtual void stripLeadingWhitespace(std::string& s) const;

    virtual bool skipLine(const std::string& s) const;


  protected:
    std::istream& _is;
    unsigned int _lineno;
    std::string _name;
    std::string _current_line;
    std::list<std::string> _lines;

  };
}




#endif
