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

// @cond PACKAGES_INTERNAL

#if !defined(LOOS_FIDLIB_HPP)
#define LOOS_FIDLIB_HPP



#include <loos.hpp>


typedef std::vector<int>                   vecInt;
typedef std::vector<uint>                 vecUint;
typedef std::vector<loos::AtomicGroup>   vecGroup;
typedef std::vector<double>             vecDouble;




// Return indices of non-zero entries in the vector (i.e. frames that are not assigned)
vecUint findFreeFrames(const vecInt& map);

// Given a set of reference structures and a trajectory, classify the trajectory
// based on which reference structure is closest to each trajectory frame
vecUint assignStructures(loos::AtomicGroup& model, loos::pTraj& traj, const vecUint& frames, const vecGroup& refs);

// Given a vector that contains indices into a trajectory, will trim off the
// end so the # of frames is an even multiple of the requested bin size (via frac)
vecUint trimFrames(const vecUint& frames, const double frac);

// Randomly partition trajectory space
// f = the fractional bin size (i.e. probability)
boost::tuple<vecGroup, vecUint> pickFiducials(loos::AtomicGroup& model, loos::pTraj& traj, const vecUint& frames, const double f);

// Find the max value in the vector
int findMaxBin(const vecInt& assignments);

// Simple histogram of bin assignments
vecUint histogramBins(const vecInt& assignments);

#endif


// @endcond TOOLS_INTERNAL
