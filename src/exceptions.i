%include <std_string.i>

%{


#include <sstream>
#include <exception>
#include <string>

#include <Atom.hpp>

%}

%include "exceptions.hpp"

%extend loos::LOOSError {
   std::string __str__() { return(self->what()); }
}

%extend loos::UnsetProperty {
   std::string __str__() { return(self->what()); }
}

%extend loos::NumericalError {
   std::string __str__() { return(self->what()); }
}

%extend loos::FileError {
   std::string __str__() { return(self->what()); }
}

%extend loos::FileOpenError {
   std::string __str__() { return(self->what()); }
}

%extend loos::FileReadError {
   std::string __str__() { return(self->what()); }
}

%extend loos::FileReadErrorWithLine {
   std::string __str__() { return(self->what()); }
}


%extend loos::FileWriteError {
   std::string __str__() { return(self->what()); }
}

%extend loos::ParseError {
   std::string __str__() { return(self->what()); }
}
