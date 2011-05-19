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


#include <sstream>
#include <exception>
#include <string>

#include <Atom.hpp>

namespace loos {


  //! Generic LOOS exception
  class LOOSError : public std::exception {
  protected:
    std::string _msg;
  public:
    explicit LOOSError() : _msg("LOOS Error") { }
    explicit LOOSError(const std::string& arg) : _msg(arg) { }
    explicit LOOSError(const Atom& a, const std::string& arg) {
      std::stringstream ss;
      ss << a << std::endl << arg;
      _msg = ss.str();
    }

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
    explicit MissingProperty(const Atom& a, const std::string& arg) : LOOSError(a, arg) { }
  };

  //! Exception caused by a blas/atlas error
  class NumericalError : public LOOSError {
  public:
    explicit NumericalError(const std::string& arg, const int info) {
      std::stringstream ss;
      ss << arg << ", info = " << info;
      _msg = ss.str();
    }
    explicit NumericalError(const std::string& arg) : LOOSError(arg) { }
  };


  //! Exception caused by inability to assign atomic numbers
  class UnknownAtomicMass : public LOOSError {
  public:
    explicit UnknownAtomicMass(const std::string& arg) : LOOSError(arg) { }
  };
};


#endif
