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

#if !defined(LOOS_HBONDDETECTOR_HPP)
#define LOOS_HBONDDETECTOR_HPP

#include <loos_defs.hpp>
#include <Coord.hpp>
#include <AtomicGroup.hpp>
#include <PeriodicBox.hpp>

namespace loos {

    //! Class for detecting hydrogen bonds
    /**
     *  This is a class intended to make it easy to determine if 
     *  3 atoms form a hydrogen bond based on distance and angle criteria.
     *  The basic criteria are set when the class is constructed (as is the 
     *  periodicity information) so that one only needs to pass the pAtoms 
     *  when using it.  The periodicity information is extracted from the 
     *  AtomicGroup passed to the constructor, so that group's shared periodic box 
     *  will be used in all distance calculations.
     *
     *  If you use the constructors with default distances and cutoffs, you get 
     *  a distance cutoff of 3.5 Ang, and an angle cutoff of 20 degrees from
     *  linear.
     *
     */
    class HBondDetector {
    public:

        HBondDetector(const double distance, const double angle, 
                      const AtomicGroup &group);
        
        HBondDetector(const AtomicGroup &group);
        
        HBondDetector();
        
        //! Method to test if this triple of atoms forms an h-bond
        /**
         *  Assumes that donor and hydrogen are in the same molecule, so 
         *  no imaging is performed when calculating their vector.  First tests
         *  if the hydrogen-acceptor distance is less than the threshold, then
         *  tests the angle to make sure it's straighter than the threshold. 
         *
         *  Actually operates on the distance**2 and the cosine of the angle, for
         *  better performance.
         */
        bool hBonded(const pAtom donor, const pAtom hydrogen, 
                     const pAtom acceptor);

    private:
        SharedPeriodicBox box;
        double cutoff_dist2;
        double cutoff_cos;

        

    };
}

#endif

