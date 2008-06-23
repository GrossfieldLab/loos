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






#if !defined(LOOS_HPP)
#define LOOS_HPP

#include <Coord.hpp>
#include <boost/shared_ptr.hpp>



typedef double greal;
typedef long gint;

typedef float dcd_real;
typedef double dcd_double;

typedef Coord<double> GCoord;
typedef boost::shared_ptr<GCoord> pGCoord;


const uint kilobytes = 1024;
const uint megabytes = kilobytes * kilobytes;
const uint gigabytes = megabytes * megabytes;



#include <utils.hpp>

#endif


