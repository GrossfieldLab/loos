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
#include <cryst.hpp>
%}

%include "cryst.hpp"



namespace loos {


    %extend UnitCell {
     std::string __repr__() {
             std::ostringstream oss;
             oss << *$self;
             size_t n = oss.str().size();
             char* buf = new char[n+1];
             strncpy(buf, oss.str().c_str(), n+1);
             return(buf);
         }

     }

 }
