/*
  KernelValue.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

*/


#if !defined(KERNELVALUE_HPP)
#define KERNELVALUE_HPP



#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include <cmath>

#include <string.h>



using namespace std;

namespace loos {

  const float FLT_THRESHOLD = 1e-10;
  
  struct Value {
    enum ValueType { NONE, STRING, INT, FLOAT };
    
    ValueType type;
    union {
      string *str;
      float flt;
      int itg;
    };

    void copy(const Value& v) {
      switch(v.type) {
      case STRING:
	str = new string(*(v.str)); break;
      case INT:
	itg = v.itg; break;
      case FLOAT:
	flt = v.flt; break;
      case NONE:
      default:
	;
      }
      type = v.type;
    }


    Value() : type(NONE) { }
    ~Value() { if (type == STRING) delete str; }

    Value(const Value& v) { copy(v); }

    const Value& operator=(const Value& v) { copy(v); return(*this); }

    Value(const string s) { setString(s); }
    Value(const float f) { setFloat(f); }
    Value(const int i) { setInt(i); }
	

    void setString(const string s) { str = new string(s); type = STRING; }
    void setFloat(const float f) { flt = f; type = FLOAT; }
    void setInt(const int i) { itg = i; type = INT; }

    string getString(void) const {
      if (type != STRING)
	throw(runtime_error("Expected a string value..."));
      return(*str);
    }

    float getFloat(void) const {
      if (type != FLOAT)
	throw(runtime_error("Expected a float value..."));
      return(flt);
    }

    int getInt(void) const {
      if (type != INT)
	throw(runtime_error("Expected an int value..."));
      return(itg);
    }

    friend ostream& operator<<(ostream& os, const Value& v) {
      switch(v.type) {
      case STRING:
	os << "<VALUE TYPE='STRING'>" << v.str << "</VALUE>";
	break;
      case FLOAT:
	os << "<VALUE TYPE='FLOAT'>" << v.flt << "</VALUE>";
	break;
      case INT:
	os << "<VALUE TYPE='INT'>" << v.itg << "</VALUE>";
	break;
      case NONE:
      default:
	os << "<VALUE TYPE='NONE'/>";
      }
      return(os);
    }

  };

  // There really oughtta be a better way to handle this...
  int compare(const Value&, const Value&);



};




#endif
