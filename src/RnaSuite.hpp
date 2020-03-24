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
#include <AtomicGroup.hpp>

using namespace std;

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

        RnaSuite(const AtomicGroup &group, const double suiteness_cutoff_);
        
        RnaSuite(const AtomicGroup &group);
        
        RnaSuite();
        
        //! Method to assign residues to backbone suites from Richardson et al.
        /**
         *  This method assigns residues to one of the 46 backbone suites
         *  defined in Richardson et al. (2008) RNA 14, 465-481. The suite of a
         *  residue is defined from delta of the previous residue to delta of
         *  the current residue.
         */
        void assignRichardsonSuites();

        //! Method to calculate backbone dihedrals for each RNA residue
        /**
         *  This method calculates the six RNA backbone dihedrals (i.e. alpha,
         *  beta, gamma, delta, epsilon, and zeta) for each residue.
         */
        void calculateBackboneDihedrals();

        //! Method to extract RNA backbone atoms from an AtomicGroup
        /**
         *  This method selects RNA backbone atoms (i.e. P, O5', C5', C4', C3',
         *  and O3') and splits them into AtomicGroups by residue id.
         */
        void extractRnaBackboneAtoms(const AtomicGroup &group);

        //! Method to return the cutoff for the suiteness score of non-outliers
        double getSuitenessCutoff() const;

        //! Method to print groups of backbone atoms for each dihedral
        void printBackboneAtoms() const;

        //! Method to set the cutoff for the suiteness score of non-outliers
        void setSuitenessCutoff(const double suiteness_cutoff_);

    private:

        std::vector<std::vector<AtomicGroup>> alpha_atoms;
        std::vector<std::vector<AtomicGroup>> beta_atoms;
        std::vector<std::vector<AtomicGroup>> gamma_atoms;
        std::vector<std::vector<AtomicGroup>> delta_atoms;
        std::vector<std::vector<AtomicGroup>> epsilon_atoms;
        std::vector<std::vector<AtomicGroup>> zeta_atoms;
        std::vector<std::vector<double>> alpha;
        std::vector<std::vector<double>> beta;
        std::vector<std::vector<double>> gamma;
        std::vector<std::vector<double>> delta;
        std::vector<std::vector<double>> epsilon;
        std::vector<std::vector<double>> zeta;
        std::vector<char> suite_name_hemi5;
        std::vector<char> suite_name_hemi3;
        std::vector<double> suiteness;
        double suiteness_cutoff;

    };
}

#endif

