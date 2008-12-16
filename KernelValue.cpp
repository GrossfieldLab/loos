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




#include <KernelValue.hpp>

namespace loos {

  namespace internal {
    


    std::string Value::getString(void) const {
      if (type != Value::STRING)
        throw(std::runtime_error("Expected a string value..."));
      return(*str);
    }
    
    float Value::getFloat(void) const {
      if (type != Value::FLOAT)
        throw(std::runtime_error("Expected a float value..."));
      return(flt);
    }

    int Value::getInt(void) const {
      if (type != Value::INT)
        throw(std::runtime_error("Expected an int value..."));
      return(itg);
    }

    //! Output in pseudo-XML
    std::ostream& operator<<(std::ostream& os, const Value& v) {
      switch(v.type) {
      case Value::STRING:
        os << "<VALUE TYPE='STRING'>" << *(v.str) << "</VALUE>";
        break;
      case Value::FLOAT:
        os << "<VALUE TYPE='FLOAT'>" << v.flt << "</VALUE>";
        break;
      case Value::INT:
        os << "<VALUE TYPE='INT'>" << v.itg << "</VALUE>";
        break;
      case Value::NONE:
      default:
        os << "<VALUE TYPE='NONE'/>";
      }
      return(os);
    }


    void Value::copy(const Value& v) {
      switch(v.type) {
      case Value::STRING:
        str = new std::string(*(v.str)); break;
      case Value::INT:
        itg = v.itg; break;
      case Value::FLOAT:
        flt = v.flt; break;
      case Value::NONE:
      default:
        ;
      }
      type = v.type;
    }
    

  
    int compare(const Value& x, const Value& y) {
      float d;
      int e;

      if (x.type != y.type)
        throw(std::runtime_error("Comparing values with different types."));

      switch(x.type) {
      case Value::STRING:
        e = (*(x.str) == *(y.str));
        if (!e) {
          if (*(x.str) < *(y.str))
            return(-1);
          else
            return(1);
        }
        return(0);

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

      case Value::NONE:
      default:
        throw(std::runtime_error("Invalid comparison"));

      }

    }
  }
}
