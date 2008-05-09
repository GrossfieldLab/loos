/*
  KernelValue.cpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

*/


#include "KernelValue.hpp"

namespace loos {

  int compare(const Value& x, const Value& y) {
    float d;
    int e;

    if (x.type != y.type)
      throw(runtime_error("Comparing values with different types."));

    switch(x.type) {
    case Value::STRING:
      return(strcmp(x.str, y.str));

    case Value::FLOAT:
      d = x.flt - y.flt;
      if (fabs(d) <= FLT_THRESHOLD)
	return(0);
      if (d < 0.0)
	return(-1);
      return(1);
      
    case Value::INT:
      e = x.itg - y.itg;
      return(e);
    }

  }


};
