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

#include <HBondDetector.hpp>

namespace loos {
    HBondDetector::HBondDetector(const double distance, const double angle, 
                                 const AtomicGroup &group) {
        
        cutoff_cos = cos((angle)*M_PI/180.0);
        cutoff_dist2 = distance * distance;
        box = group.sharedPeriodicBox();
    }

    HBondDetector::HBondDetector(const AtomicGroup &group) {
        cutoff_cos = cos((20.)*M_PI/180.);
        cutoff_dist2 = 3.5 * 3.5; 
        box = group.sharedPeriodicBox();
    }
        

    HBondDetector::HBondDetector() {
        cutoff_cos = cos((20.)*M_PI/180.);
        cutoff_dist2 = 3.5 * 3.5; 
        box = SharedPeriodicBox();
    }
        

    bool HBondDetector::hBonded(const pAtom donor, const pAtom hydrogen,
                                const pAtom acceptor) {
        // Check distance between hydrogen and acceptor
        double d2;
        if (box.isPeriodic()) {
            d2 = hydrogen->coords().distance2(acceptor->coords(), box.box());
        }
        else {
            d2 = hydrogen->coords().distance2(acceptor->coords());
        }

        if (d2 > cutoff_dist2) {
            return false;
        }


        // If the distance test passes, try the angle test
        // We assume the donor and hydrogen are in the same periodic image
        // Return true if the angle is greater than the threshold (meaning
        // the cosine is less than the threshold)
        GCoord d_to_h =  hydrogen->coords() - donor->coords();
        GCoord h_to_a =  acceptor->coords() - hydrogen->coords();
        if (box.isPeriodic()) {
            h_to_a.reimage(box.box());
        }

        double cosine = (d_to_h * h_to_a)/(d_to_h.length() * sqrt(d2));
        return (cosine > cutoff_cos);
        }

}
