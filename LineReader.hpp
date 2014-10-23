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


  //! Class for reading line-by-line from a file while tracking line numbers and stripping comments

  class LineReader {
  public:
    LineReader() : _is(0), _lineno(1) { }
    LineReader(std::istream& is) : _is(&is), _lineno(1) { }
    LineReader(std::istream& is, const std::string& name) : _is(&is), _lineno(1), _name(name) { }
    
    virtual ~LineReader() { }

    //! Access the internal stream pointer
    virtual std::istream& stream() const;

    //! Set the internal stream pointer
    virtual void stream(std::istream& is);
    
    //! Access the name associated with the internal stream
    virtual std::string name() const;

    //! Set the name associated with the internal stream
    virtual void name(const std::string& name);

    //! Get the next line from the file, returning true if successful
    virtual bool getNext();

    //! Put a line back onto the file (virtually)
    virtual void push_back(const std::string& s);

    //! The currently read line
    virtual std::string line() const;
    
    //! The current line number into the file
    virtual uint lineNumber() const;

  protected:


    /*
      The following member functions are the ones that will most
      likely need to be changed in derived classes to alter the
      functionality of the reader...
    */

    // Verifies the state of the stream (throws if error)
    virtual void checkState() const throw(FileParseError);

    // Remove comments from the string
    virtual void stripComment(std::string& s) const;

    // Remove leading whitespace (space and tabs)
    virtual void stripLeadingWhitespace(std::string& s) const;

    // True means to skip the current line (e.g. the line is empty)
    virtual bool skipLine(const std::string& s) const;


  protected:
    std::istream* _is;
    unsigned int _lineno;
    std::string _name;
    std::string _current_line;
    std::list<std::string> _lines;

  };
}




#endif
