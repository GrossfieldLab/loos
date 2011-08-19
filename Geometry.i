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



%header %{
#include <Geometry.hpp>
%}

namespace loos {
  //! Namespace for math and math-related things in loos
  namespace Math {

    const double DEGREES = 180 / M_PI ;
    
    //! Compute the angle in degrees assuming the middle is the vertex
    greal angle(const GCoord&, const GCoord&, const GCoord&, const GCoord* =NULL);
    
    //! Compute the angle in degrees assuming the middle is the vertex
    greal angle(const pAtom& a, const pAtom& b, const pAtom& c, const GCoord * =NULL);
    
    //! Compute the torsion in degrees 
    greal torsion(const GCoord& a, const GCoord& b, const GCoord& c, 
                  const GCoord& d, const GCoord* =NULL);
    
    //! Compute the torsion in degrees 
    greal torsion(const pAtom& a, const pAtom& b, const pAtom& c, 
                  const pAtom& d, const GCoord* =NULL);

  }
}
