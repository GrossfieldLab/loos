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


%include <std_string.i>

%header %{
#include <XForm.hpp>

  typedef double                   greal;
  typedef loos::Matrix44<double>   GMatrix;
  typedef loos::Coord<double>      GCoord;

%}


%include "XForm.hpp"

namespace loos {

   %rename(XFormVector)    std::vector<XForm>;
};



%extend loos::XForm {


  std::string __repr__() {
    static char buf[1024];

    GMatrix M = $self->current();
    std::ostringstream oss;
    oss << "[ ";
    for (uint j=0; j<4; ++j) {
      oss << "[";
      for (uint i=0; i<4; ++i)
	oss << M(j, i) << (i < 3 ? ", " : " ], ");
    }
    oss << "]";
    return(oss.str());
  }


  

};
