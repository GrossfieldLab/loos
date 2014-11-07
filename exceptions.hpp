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


#if !defined(LOOS_EXCEPTIONS_HPP)
#define LOOS_EXCEPTIONS_HPP


#include <sstream>
#include <exception>
#include <string>
#include <loos_defs.hpp>

namespace loos {


  // Forward decl's
  class Atom;
  std::ostream& operator<<(std::ostream&, const Atom&);

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

  //! Exception in options
  class OptionsError : public LOOSError 
  {
  public:
    explicit OptionsError(const std::string& arg) : LOOSError(arg) { }
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

  //! Exception caused by insufficient atom properties.. (DEPRECATED)
  class MissingProperty : public LOOSError {
  public:
    explicit MissingProperty(const std::string& arg) : LOOSError(arg) { }
    explicit MissingProperty(const Atom& a, const std::string& arg) : LOOSError(a, arg) { }
  };

  class UnsetProperty : public LOOSError {
  public:
    UnsetProperty() : LOOSError("Attempting to access an unset atom property") {}
    UnsetProperty(const std::string& p) : LOOSError(p) {}
    UnsetProperty(const Atom& a, const std::string& p) : LOOSError(a, p) {}
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

  class FileError : public LOOSError {
  protected:
    std::string _operation;
    std::string _filename;
  public:
    FileError(const std::string& op) : LOOSError("Error while " + op), _operation(op) {}
    FileError(const std::string& op, const std::string& fname)
      : LOOSError("Error while " + op + " from " + fname),
        _operation(op), _filename(fname)
    {}

    FileError(const std::string& op,
              const std::string& fname,
              const std::string& msg)
      : LOOSError("Error while " + op + " from " + fname + msg),
        _operation(op),
        _filename(fname)
    {}

    std::string operation() const throw() { return(_operation); }
    std::string filename() const throw() { return(_filename); }


    ~FileError() throw() {}
  };

  class FileOpenError : public FileError {
  public:
    FileOpenError() : FileError("opening") { }
    FileOpenError(const std::string& fname) : FileError("opening", fname) {}
    FileOpenError(const std::string& fname, const std::string& msg) : FileError("opening", fname, '\n' + msg) {}
  };


  class FileReadError : public FileError {
  public:
    FileReadError() : FileError("reading from") { }
    FileReadError(const std::string& fname) : FileError("reading from", fname) {}
    FileReadError(const std::string& fname, const std::string& msg) : FileError("reading from", fname, '\n' + msg) {}
  };


  class FileReadErrorWithCode : public FileReadError {
  protected:
    uint _error_code;
    std::string _msg;
  public:
    FileReadErrorWithCode(const uint n)
      : FileReadError("reading"), _error_code(n)
    { init(); }

    FileReadErrorWithCode(const std::string& fname, const uint n)
      : FileReadError("reading", fname), _error_code(n)
    { init(); }

    FileReadErrorWithCode(const std::string& fname, const std::string& msg, const uint n)
      : FileReadError("reading", fname), _error_code(n), _msg(msg)
    { init(); }

    uint errorCode() const { return(_error_code); }
    
    ~FileReadErrorWithCode() throw() {}

  private:
    void init() {
      std::ostringstream oss;

      oss << FileReadError::_msg << ", error code " << _error_code << std::endl << _msg;
      FileReadError::_msg = oss.str();
    }

  };


  class FileReadErrorWithLine : public FileReadError {
  protected:
    uint _lineno;
    std::string _msg;
  public:
    FileReadErrorWithLine(const uint ln)
      : FileReadError("reading"), _lineno(ln)
    { init(); }

    FileReadErrorWithLine(const std::string& fname, const uint ln)
      : FileReadError("reading ", fname), _lineno(ln)
    { init(); }
    

    FileReadErrorWithLine(const std::string& fname, const std::string& msg, const uint ln)
      : FileReadError("reading ", fname), _lineno(ln), _msg(msg)
    { init(); }

    uint lineNumber() const throw() { return(_lineno); }

    ~FileReadErrorWithLine() throw() {}

  private:
    void init() {
      std::ostringstream oss;
      
      oss << FileReadError::_msg << " at line " << _lineno << std::endl << _msg;
      FileReadError::_msg = oss.str();
    }

  };




  class FileWriteError : public FileError {
  public:
    FileWriteError() : FileError("writing to") { }
    FileWriteError(const std::string& fname) : FileError("writing to", fname) {}
    FileWriteError(const std::string& fname, const std::string& msg) : FileError("writing to", fname, '\n' + msg) {}
  };
  


  //! Exception for writing trajectories
  class TrajectoryWriteError : public LOOSError {
  public:
    TrajectoryWriteError() : LOOSError("Error while writing trajectory") {}
    TrajectoryWriteError(const std::string& arg) : LOOSError(arg) {}
  };

  //! Exception for Reading trajectories
  class TrajectoryReadError : public LOOSError {
  public:
    TrajectoryReadError() : LOOSError("Error while reading from trajectory") {}
    TrajectoryReadError(const std::string& arg) : LOOSError(arg) {}
  };


  class EndOfFile : public LOOSError {
  public:
    EndOfFile() : LOOSError("Attempting to read past end of file") {}
    EndOfFile(const std::string& arg) : LOOSError(arg) {}
  };

  //! Exceptions for reading Amber NetCDF files
  struct AmberNetcdfError : public LOOSError {
    explicit AmberNetcdfError(const std::string& msg) : LOOSError(msg) { }
    explicit AmberNetcdfError(const std::string& msg, const int retval) {
      std::stringstream ss;
      ss << msg << " with error #" << retval;
      _msg = ss.str();
    }

  };

  struct AmberNetcdfOpenError : public AmberNetcdfError {
    explicit AmberNetcdfOpenError() : AmberNetcdfError("Error opening Amber NetCDF file") { }
  };

  struct AmberNetcdfTypeError : public AmberNetcdfError {
    explicit AmberNetcdfTypeError(const std::string msg) : AmberNetcdfError(msg) { }
  };
  

  struct BadConnectivityError : public LOOSError {
    explicit BadConnectivityError(const std::string& msg) : LOOSError(msg) {}
  };



  class FileParseError : public LOOSError {
  public:
    explicit FileParseError(const std::string& arg, const uint lineno) {
      std::stringstream ss;
      ss << arg << " at line " << lineno;
      _msg = ss.str();
    }
  };

};

#endif
