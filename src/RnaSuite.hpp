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

#if !defined(LOOS_RNASUITE_HPP)
#define LOOS_RNASUITE_HPP

#include <loos_defs.hpp>
#include <Coord.hpp>
#include <AtomicGroup.hpp>
#include <PeriodicBox.hpp>

namespace loos {

    //! Class for assigning backbone suites to an RNA
    /**
     *  This class acts on an AtomicGroup and assigns backbone suites (as 
     *  defined in Richardson et al. (2008) RNA 14, 465-481) to any RNA residues
     *  present. It also calculates the "suiteness" score that describes how 
     *  well the residue fits into its assigned suite.
     */
    class RnaSuite {
    public:

        RnaSuite(const AtomicGroup &group, const double &suiteness_tolerance);
        
        RnaSuite(const AtomicGroup &group);
        
        RnaSuite();
        
        //! Method to extract RNA backbone atoms from an AtomicGroup
        /**
         *  This method selects RNA backbone atoms (i.e. P, O5', C5', C4', C3', 
         *  and O3') and splits them into AtomicGroups by residue id.
         */
        void extractRnaBackboneAtoms();

        //! Method to calculate backbone dihedrals for each RNA residue
        /**
         *  This methods calculates the six RNA backbone dihedrals (i.e. alpha, 
         *  beta, gamma, delta, epsilon, and zeta) for each residue.
         */
        void calculateBackboneDihedrals();

    private:
        vector<vector<AtomicGroup>> alpha_atoms;
        vector<vector<AtomicGroup>> beta_atoms;
        vector<vector<AtomicGroup>> gamma_atoms;
        vector<vector<AtomicGroup>> delta_atoms;
        vector<vector<AtomicGroup>> epsilon_atoms;
        vector<vector<AtomicGroup>> zeta_atoms;
        double suiteness_tolerance;

    };
}

#endif

