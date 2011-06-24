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




#include "water-lib.hpp"


namespace loos {
  namespace DensityTools {

    std::vector<GCoord> getBounds(pTraj& traj, AtomicGroup& g, const std::vector<uint>& indices) {
      double d = std::numeric_limits<double>::max();
      GCoord min(d, d, d);
      GCoord max(-d, -d, -d);
      
      for (std::vector<uint>::const_iterator i = indices.begin(); i != indices.end(); ++i) {
        traj->readFrame(*i);
        traj->updateGroupCoords(g);
        std::vector<GCoord> bdd = g.boundingBox();
        for (int j=0; j<3; ++j) {
          if (bdd[0][j] < min[j])
            min[j] = bdd[0][j];
          if (bdd[1][j] > max[j])
            max[j] = bdd[1][j];
        }
      }
      
      std::vector<GCoord> bdd;
      bdd.push_back(min);
      bdd.push_back(max);
      return(bdd);
    }
 

  };

};
