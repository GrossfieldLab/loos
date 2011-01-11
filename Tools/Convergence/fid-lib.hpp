/*
  Fiducial library
*/


/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010, Tod D. Romo
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



#if !defined(FIDLIB_HPP)
#define FIDLIB_HPP



#include <loos.hpp>


typedef std::vector<int>                   vecInt;
typedef std::vector<uint>                 vecUint;
typedef std::vector<loos::AtomicGroup>   vecGroup;
typedef std::vector<double>             vecDouble;





vecUint findFreeFrames(const vecInt& map);

vecUint assignStructures(loos::AtomicGroup& model, loos::pTraj& traj, const vecUint& frames, const vecGroup& refs);
vecUint trimFrames(const vecUint& frames, const double frac);
boost::tuple<vecGroup, vecUint> pickFiducials(loos::AtomicGroup& model, loos::pTraj& traj, const vecUint& frames, const double f);

int findMaxBin(const vecInt& assignments);
vecUint histogramBins(const vecInt& assignments);



#endif
