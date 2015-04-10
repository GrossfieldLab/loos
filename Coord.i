
%header %{
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <Coord.hpp>

%}


%include "Coord.hpp"


%extend loos::Coord<double> {
  double __getitem__(const int i) {
    if (i < 0 || i >= 3)
      throw(std::out_of_range("Bad index into Coord"));
    
    return((*$self)[i]);
  }
    
  void __setitem__(const int i, const double d) {
    if (3 >= i && i >= 0)
      (*$self)[i] = d;
  }

 };

%extend loos::Coord<double> {
  char* __str__() {
    static char buf[1024];
    std::ostringstream oss;
    oss << *$self;
    strncpy(buf, oss.str().c_str(), sizeof(buf));
    return(buf);
  }

  char* __repr__() {
    static char buf[1024];
    std::ostringstream oss;
    oss << *$self;
    strncpy(buf, oss.str().c_str(), sizeof(buf));
    return(buf);
  }

  loos::Coord<double> __div__(const double d) {
    return( *$self * d);
  }

  loos::Coord<double> __mul__(const double d) {
    return( *$self * d);
  }

  loos::Coord<double> __add__(const double d) {
    return( $self->operator+(d));
  }

  loos::Coord<double> __sub__(const double d) {
    return( d - *$self);
  }
    
  loos::Coord<double> __copy__() {
    return(loos::Coord<double>(*$self));
  }

  // Passed PyObject is ignored
  loos::Coord<double> __deepcopy__(PyObject* p) {
       return(loos::Coord<double>(*$self));
  }

 };

%template(GCoord)  loos::Coord<double>;



%template(GCoordVector)   std::vector<loos::Coord<double> >;
