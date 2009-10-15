/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo, Alan Grossfield
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


#if !defined(EXCEPTIONS_HPP)
#define EXCEPTIONS_HPP



#include <exception>
#include <string>



namespace loos {


  //! Generic LOOS exception
  class LOOSError : public std::exception {
    std::string _msg;
  public:
    explicit LOOSError(const std::string& arg) : _msg(arg) { }
    virtual ~LOOSError() throw() {};
    virtual const char* what(void) const throw() { return(_msg.c_str()); }
  };

  
  //! Exception when parsing input data
  class ParseError : public LOOSError {
  public:
    explicit ParseError(const std::string& arg) : LOOSError(arg) { }
  };

  //! Exception cause by some operation failing (ie no atoms selected)
  class NullResult : public LOOSError {
  public:
    explicit NullResult(const std::string& arg) : LOOSError(arg) { }
  };

  //! Exception caused by insufficient atom properties..
  class MissingProperty : public LOOSError {
  public:
    explicit MissingProperty(const std::string& arg) : LOOSError(arg) { }
  };



};


#endif
