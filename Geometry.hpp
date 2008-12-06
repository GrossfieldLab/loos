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




#if !defined(GEOM_HPP)
#define GEOM_HPP

#include "loos_defs.hpp"

#include <cmath>
#include "Atom.hpp"

namespace loos {
  namespace Math {

    const double DEGREES = 180 / M_PI ;
    
    //! Compute the angle in degrees assuming the middle is the vertex
    greal angle(const GCoord &a, const GCoord &b, const GCoord &c);
    
    //! Compute the angle in degrees assuming the middle is the vertex
    greal angle(const pAtom a, const pAtom b, const pAtom c);
    
    //! Compute the torsion in degrees 
    greal torsion(const GCoord &a, const GCoord &b, const GCoord &c, 
                  const GCoord &d);
    
    greal torsion(const pAtom a, const pAtom b, const pAtom c, 
                  const pAtom d);
  }
}

#endif
