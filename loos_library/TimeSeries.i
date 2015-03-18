/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Alan Grossfield
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


%header %{
#include <TimeSeries.hpp>
%}

%include "TimeSeries.hpp"

namespace loos {

  %extend TimeSeries<double> {

    double __getitem__(const int i) {
      if (i < 0 || static_cast<uint>(i) >= $self->size())
	throw(std::out_of_range("Index into TimeSeries is out of bounds"));
      
      return((*$self)[i]);
    }
    
    void __setitem__(const int i, const double d) {
      (*$self)[i] = d;
    }
  };

  %rename(__add__)  loos::TimeSeries<double>::operator+;
  %rename(__sub__) loos::TimeSeries<double>::operator-;
  %rename(__mul__) loos::TimeSeries<double>::operator*;
  %rename(__div__) loos::TimeSeries<double>::operator/;


  %template(TimeSeriesDbl) loos::TimeSeries<double>;
  



}
