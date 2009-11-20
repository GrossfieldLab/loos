/*
  enm-lib

  (c) 2009 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry


  Common code for the ENM suite

*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009 Tod D. Romo
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



#if !defined(ENMLIB_HPP)
#define ENMLIB_HPP

#include <loos_defs.hpp>
#include <GCoord.hpp>
#include <AtomicGroup.hpp>
#include <Selectors.hpp>

#include <vector>


std::vector<loos::AtomicGroup> sideChainCentroids(const loos::AtomicGroup& grp, int maxid = 0, int maxresid = 0,
                                            const std::string& name = "SID", const std::string& resname = "SID", const std::string& segid = "SIDE");





#endif
