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
    explicit LOOSError(const std::string& fname, const Atom& a, const std::string& arg) {
        std::stringstream ss;
        ss << "In file: " + fname << std::endl << a << std::endl << arg;
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


  //! Exception when trying to use an unset Atom property
  class UnsetProperty : public LOOSError {
  public:
    UnsetProperty() : LOOSError("Attempting to access an unset atom property") {}
    UnsetProperty(const std::string& p) : LOOSError(p) {}
    UnsetProperty(const Atom& a, const std::string& p) : LOOSError(a, p) {}
  };


  //! Exception indicating internal XDR error
  /**
   * This most likely means that the word-size of your data
   * exceeds what the LOOS XDR lib was built to handle
   */
  struct XDRDataSizeError : public LOOSError {
    XDRDataSizeError() : LOOSError("XDR data size error") {}
    XDRDataSizeError(const std::string& s) : LOOSError(s) {}
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


  //! Errors related to File I/O
  /**
   * Most file I/O exceptions derive from this class.
   *
   * Note that most I/O classes currently do not guarantee 
   * that they will be in a safe/usable state once an 
   * exception occurs
   */
  class FileError : public LOOSError {
  protected:
    std::string _operation;
    std::string _filename;
    int _errcode;

  public:
    FileError(const std::string& op) : LOOSError("Error while " + op), _operation(op) {}

    FileError(const std::string& op, const std::string& fname)
      : LOOSError("Error while " + op + " " + fname),
        _operation(op), _filename(fname)
    {}

    FileError(const std::string& op,
              const std::string& fname,
              const std::string& msg)
      : LOOSError("Error while " + op + " " + fname + msg),
        _operation(op),
        _filename(fname)
    {}

    FileError(const std::string& op,
              const std::string& fname,
              const std::string& msg,
              const int err)
      : LOOSError("Error while " + op + " " + fname + msg),
        _operation(op),
        _filename(fname),
        _errcode(err)
    {}



    //! What operation was involved (e.g. reading, writing. etc)
    std::string operation() const throw() { return(_operation); }

    //! File that had the problem (or "stream" if not a file)
    std::string filename() const throw() { return(_filename); }

    //! The error code that may have been generated
    int errorCode() const { return(_errcode); }

    //! Sets the error code
    void errorCode(const int i) { _errcode = i; }

    ~FileError() throw() {}
  };

  //! Error while opening a file
  /**
   * This mostly represents issues with actually opening the file, such as
   * a bad filename, permissions, as well as scanning the file as part of
   * instantiation.
   *
   * Some classes will automatically read either the whole file or, in the
   * case of trajectories, the first frame.  This can result in a FileReadError
   * exception being thrown, despite being part of a logical "open" operation.
   *
   */
  class FileOpenError : public FileError {
  public:
    FileOpenError() : FileError("opening") { }
    FileOpenError(const std::string& fname) : FileError("opening", fname) {}
    FileOpenError(const std::string& fname, const std::string& msg) : FileError("opening", fname, '\n' + msg) {}
    FileOpenError(const std::string& fname, const std::string& msg, const int err) : FileError("opening", fname, '\n' + msg, err) {}
  };


  //! Errors that occur while reading a file
  class FileReadError : public FileError {
  public:
    FileReadError() : FileError("reading from") { }
    FileReadError(const std::string& fname) : FileError("reading from", fname) {}
    FileReadError(const std::string& fname, const std::string& msg) : FileError("reading from", fname, '\n' + msg) {}
    FileReadError(const std::string& fname, const std::string& msg, const int err) : FileError("reading from", fname, '\n' + msg, err) {}
  };


  //! Errors that occur while reading a text file (where lines are tracked)
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

    //! The line number that caused the problem
    uint lineNumber() const throw() { return(_lineno); }

    ~FileReadErrorWithLine() throw() {}

  private:
    void init() {
      std::ostringstream oss;
      
      oss << FileReadError::_msg << " at line " << _lineno << std::endl << _msg;
      FileReadError::_msg = oss.str();
    }

  };

  
  //! Errors while writing to files
  class FileWriteError : public FileError {
  public:
    FileWriteError() : FileError("writing to") { }
    FileWriteError(const std::string& fname) : FileError("writing to", fname) {}
    FileWriteError(const std::string& fname, const std::string& msg) : FileError("writing to", fname, '\n' + msg) {}
  };
  

  //! Errors related to trajectory reading and writing
  /**
   * Most trajectory exceptions derive from this class.
   */
  class TrajectoryError : public LOOSError {
  protected:
      std::string _operation;
      std::string _filename;
      int _errcode;

  public:
      TrajectoryError(const std::string& op) : LOOSError("Error while " + op), _operation(op) {}

      TrajectoryError(const std::string& op, const std::string& fname)
          : LOOSError("Error while " + op + ", " + fname),
          _operation(op), _filename(fname)
      {}

      TrajectoryError(const std::string& op,
          const std::string& fname,
          const std::string& msg)
          : LOOSError("Error while " + op + ", " + fname + "\n" + msg),
          _operation(op),
          _filename(fname)
      {}

      TrajectoryError(const std::string& op,
          const std::string& fname,
          const std::string& msg,
          const int err)
          : LOOSError("Error while " + op + ", " + fname + "\n" + msg),
          _operation(op),
          _filename(fname),
          _errcode(err)
      {}

      
      //! What operation was involved (e.g. reading, writing. etc)
      std::string operation() const throw() { return(_operation); }

      //! File that had the problem (or "stream" if not a file)
      std::string filename() const throw() { return(_filename); }

      //! The error code that may have been generated
      int errorCode() const { return(_errcode); }

      //! Sets the error code
      void errorCode(const int i) { _errcode = i; }

      ~TrajectoryError() throw() {}
  };



};

#endif
