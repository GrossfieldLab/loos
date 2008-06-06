/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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

  const float FLT_THRESHOLD = 1e-10;     // Threshold for floating equality

  //! Value class for the LOOS Kernel (virtual machine)
  /**
   * These are the "values" that are on the data
   * stack for the Kernel...
  */

  struct Value {
    //! Type of data this Value contains...
    enum ValueType { NONE, STRING, INT, FLOAT };  // What type of data
						  // the union contains...
    
    ValueType type;
    union {
      string *str;
      float flt;
      int itg;
    };

    //! Make a copy (clone) of a Value.
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

    //! Retrieve data, throwing an error if the Value is of the incorrect type.
    string getString(void) const {
      if (type != STRING)
	throw(runtime_error("Expected a string value..."));
      return(*str);
    }

    //! Retrieve data, throwing an error if the Value is of the incorrect type.
    float getFloat(void) const {
      if (type != FLOAT)
	throw(runtime_error("Expected a float value..."));
      return(flt);
    }

    //! Retrieve data, throwing an error if the Value is of the incorrect type.
    int getInt(void) const {
      if (type != INT)
	throw(runtime_error("Expected an int value..."));
      return(itg);
    }

    //! Output in pseudo-XML
    friend ostream& operator<<(ostream& os, const Value& v) {
      switch(v.type) {
      case STRING:
	os << "<VALUE TYPE='STRING'>" << *(v.str) << "</VALUE>";
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


  //! Compare two Value objects, depending on their types.
  /**
   * Returns -1 if x < y
   *          0 if x = y
   *         1 if x > y
   *
   * For strings, this is based on the lexical value...
   * For floats, the equality is determined by the FLT_THRESHOLD
   * constant...
   *
   * Comparison of non-like types is an error...
   */
  int compare(const Value&, const Value&);



};




#endif
