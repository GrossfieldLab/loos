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


#include <loos.hpp>
#include <vector>

#include <boost/tuple/tuple.hpp>

namespace loos {
  //! Compute the average structure of a set of AtomicGroup objects
  /**Is xform-aware */
  AtomicGroup averageStructure(const vector<AtomicGroup>& ensemble);

  //! Compute an iterative superposition (a la Alan)
  /**Is xform-aware */
  boost::tuple<vector<XForm>, greal, int> iterativeAlignment(vector<AtomicGroup>& ensemble, greal threshold, int maxiter=1000);
};



#endif
