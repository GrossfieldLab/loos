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


#include <exceptions.hpp>
#include <LineReader.hpp>

namespace loos {


   bool LineReader::getNext() {
    if (! _lines.empty()) {
      _current_line = _lines.back();
      _lines.pop_back();
    } else {
      while (getline(_is, _current_line).good() ) {
        ++_lineno;
        stripComment(_current_line);
        stripLeadingWhitespace(_current_line);
        if (! skipLine(_current_line) )
          break;
      }
    }
    
    checkState();
    return( _is.good() );
  }

   void LineReader::push_back(const std::string& s) {
    _lines.push_back(s);
  }


   bool LineReader::eof() const { return(_is.eof()); }
   bool LineReader::fail() const { return(_is.fail()); }
   bool LineReader::good() const { return(_is.good()); }
  
   std::string LineReader::line() const { return(_current_line); }
  
   unsigned int LineReader::lineNumber() const { return(_lineno); }

   void LineReader::checkState() const {
     if (_is.fail()) {
       if (_name.empty())
         throw(FileParseError("Error while reading from " + _name, _lineno));
       else
         throw(FileParseError("Error while reading file", _lineno));
     }
   }

  void LineReader::stripComment(std::string& s) const {
    std::string::size_type i = s.find('#');
    if (i != std::string::npos)
      s.erase(i, s.length() - i);
  }
  
  void LineReader::stripLeadingWhitespace(std::string& s) const {
    std::string::size_type i = s.find_first_not_of(" \t");
    if (i > 0)
      s.erase(0, i);
  }
  
  bool LineReader::skipLine(const std::string& s) const {
    return(s.empty());
  }
  

}
