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



#if !defined(ENSEMBLES_HPP)
#define ENSEMBLES_HPP

#include <vector>
#include <XForm.hpp>

#include <loos_defs.hpp>

namespace loos {
  //! Compute the average structure of a set of AtomicGroup objects
  AtomicGroup averageStructure(const vector<AtomicGroup>& ensemble);

  //! Compute the average structure by reading through a trajectory
  AtomicGroup averageStructure(const AtomicGroup&, const vector<XForm>&, pTraj);

  //! Compute an iterative superposition (a la Alan)
  boost::tuple<vector<XForm>, greal, int> iterativeAlignment(vector<AtomicGroup>& ensemble, greal threshold, int maxiter=1000);

  //! Compute an iterative superposition by reading in frames from the Trajectory.
  /*!
   * This function will internally cache an AtomicGroup copy for each
   * frame of the trajectory.  This could chew up a lot of memory, but
   * we make the assumption that you will usually be aligning against
   * a fairly small subset of each frame...
   */
  boost::tuple<vector<XForm>, greal, int> iterativeAlignment(const AtomicGroup& g, pTraj, greal threshold, int maxiter=1000);
};



#endif
