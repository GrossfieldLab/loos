// -------------------------------------------------
// Water (density) Library
// -------------------------------------------------

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009 Tod D. Romo, Alan Grossfield
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




#if !defined(WATER_LIB_HPP)
#define WATER_LIB_HPP


#include <loos.hpp>


namespace loos {
  namespace DensityTools {

    //! Get the max bounding box for a group over the trajectory
    std::vector<GCoord> getBounds(pTraj& traj, AtomicGroup& group, const std::vector<uint>& frames);

  };

};





#endif
