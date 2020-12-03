%{


#include <sstream>
#include <exception>
#include <string>

#include <Atom.hpp>

%}

%include "exceptions.hpp"

%extend loos::LOOSError {
  const char* __str__() { std::string s = self->what(); return(s.c_str()); }
}

%extend loos::UnsetProperty {
  const char* __str__() { std::string s = self->what(); return(s.c_str()); }
}

%extend loos::NumericalError {
  const char* __str__() { std::string s = self->what(); return(s.c_str()); }
}

%extend loos::FileError {
  const char* __str__() { std::string s = self->what(); return(s.c_str()); }
}

%extend loos::FileOpenError {
  const char* __str__() { std::string s = self->what(); return(s.c_str()); }
}

%extend loos::FileReadError {
  const char* __str__() { std::string s = self->what(); return(s.c_str()); }
}

%extend loos::FileReadErrorWithLine {
  const char* __str__() { std::string s = self->what(); return(s.c_str()); }
}


%extend loos::FileWriteError {
  const char* __str__() { std::string s = self->what(); return(s.c_str()); }
}

%extend loos::ParseError {
  const char* __str__() { std::string s = self->what(); return(s.c_str()); }
}
