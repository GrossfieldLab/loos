/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2014, Tod D. Romo, Alan Grossfield
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


/*
  NOTE:  This code is based on the xdrfile library authored by:
    David van der Spoel <spoel@gromacs.org>
    Erik Lindahl <lindahl@gromacs.org>
  and covered by the GLPL-v3 license
*/
#include <xdr.hpp>



template<> uint loos::internal::XDRReader::read<double>(double* p) 
{
  double result;
  stream->read(reinterpret_cast<char*>(&result), sizeof(double));
  if (need_to_swab)
    result = swab(result);
	    
  *p = result;
  return(!stream->fail());
}


template<> uint loos::internal::XDRWriter::write<double>(const double* p) {

  unsigned long u;
  double* up = reinterpret_cast<double*>(&u);
  *up = *p;

  if (need_to_swab)
    u = swab(u);

    
  stream->write(reinterpret_cast<char*>(&u), sizeof(double));

  return(!stream->fail());
}
