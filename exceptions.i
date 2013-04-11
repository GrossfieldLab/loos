%{
    

#include <sstream>
#include <exception>
#include <string>

#include <Atom.hpp>

%}

namespace loos {


  //! Generic LOOS exception
  class LOOSError : public std::exception {
  protected:
    std::string _msg;
  public:
      explicit LOOSError();
      explicit LOOSError(const std::string& arg);
      explicit LOOSError(const Atom& a, const std::string& arg);
      virtual ~LOOSError() throw();
      virtual const char* what(void) const throw();
  };

  
  //! Exception when parsing input data
  class ParseError : public LOOSError {
  public:
      explicit ParseError(const std::string& arg);
      
  };

  class FileParseError : public LOOSError {
  public:
      explicit FileParseError(const std::string& arg, const uint lineno);

  };

  //! Exception cause by some operation failing (ie no atoms selected)
  class NullResult : public LOOSError {
  public:
      explicit NullResult(const std::string& arg);
      
  };

  //! Exception caused by insufficient atom properties..
  class MissingProperty : public LOOSError {
  public:
      explicit MissingProperty(const std::string& arg);
      explicit MissingProperty(const Atom& a, const std::string& arg);
  };

  //! Exception caused by a blas/atlas error
  class NumericalError : public LOOSError {
  public:
      explicit NumericalError(const std::string& arg, const int info);
      explicit NumericalError(const std::string& arg);
      
  };


  //! Exception caused by inability to assign atomic numbers
  class UnknownAtomicMass : public LOOSError {
  public:
      explicit UnknownAtomicMass(const std::string& arg);
  };
};
