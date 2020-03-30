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
#include <Geometry.hpp>

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
        
        //! Method to assign residues to a delta(i-1), delta, gamma index
        size_t assignDDGIndex(double dihedral, vector<double> &min,
            vector<double> &max, uint increment, uint &ddg_index);

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

        //! Calculate a dihedral in deg from 4 atoms in the range [0, 360]
        double calculateDihedral(const AtomicGroup &group);

        //! Method to check the size of a vector of continuous groups
        void checkContinuousGroupSize(
            const vector<vector<AtomicGroup>> &group_vector,
            const size_t target_size, const string dihedral_name) const;

        //! Method to check the size of a vector of residues
        void checkResidueSize(const vector<AtomicGroup> &residue_vector,
            const size_t target_size, const string dihedral_name,
            const size_t group_index) const;

        //! Method to define suites used for assignment from an existing scheme
        void defineSuites(const string suite_definition);

        //! Method to define suites used for assignment from a file
        void defineSuitesFromFile(const string filename);

        //! Method to define suites used for assignment from suitename
        void defineSuitesFromSuitename();

        //! Method to extract RNA backbone atoms from an AtomicGroup
        /**
         *  This method selects RNA backbone atoms (i.e. P, O5', C5', C4', C3',
         *  and O3') and splits them into AtomicGroups by residue id.
         */
        void extractRnaBackboneAtoms(const AtomicGroup &group);

        //! Method to return the cutoff for the suiteness score of non-outliers
        double getSuitenessCutoff() const;

        //! Calculate a scaled hyperellipsoid distance between two points
        double hyperellipsoidDist(vector<double> &dihedrals,
            vector<double> &reference, uint first_index, uint last_index);

        //! Method to print groups of backbone atoms for each dihedral
        void printBackboneAtoms() const;

        //! Method to print backbone dihedrals for each residue
        void printBackboneDihedrals() const;

        //! Method to print reference suite names and mean dihedrals
        void printReferenceSuites() const;

        //! Method to print suite names, suiteness scores, and dihedrals
        void printSuites() const;

        //! Method to set the cutoff for the suiteness score of non-outliers
        void setSuitenessCutoff(const double suiteness_cutoff_);

    private:

        // Reference suites used for assignment
        vector<vector<vector<double>>> reference_suite_dihedrals;
        vector<vector<string>> reference_suite_names;
        vector<string> reference_suite_ddgs;

        // Widths used to scale each dihedral dimension
        vector<double> dihedral_width;

        // Satellite widths used to scale overlapping clusters
        vector<double> satellite_width;

        // Boundaries for allowed regions of delta(i-1), delta, and gamma
        vector<double> delta_min;
        vector<double> delta_max;
        vector<double> gamma_min;
        vector<double> gamma_max;

        // Boundaries used to filter suites based on epsilon, zeta, alpha, beta
        vector<double> filter_min;
        vector<double> filter_max;

        // Vector of continuous groups, composed of vectors of AtomicGroups
        // for each residue within a continuous group
        vector<vector<AtomicGroup>> alpha_atoms;
        vector<vector<AtomicGroup>> beta_atoms;
        vector<vector<AtomicGroup>> gamma_atoms;
        vector<vector<AtomicGroup>> delta_atoms;
        vector<vector<AtomicGroup>> epsilon_atoms;
        vector<vector<AtomicGroup>> zeta_atoms;

        // Vector of vectors of backbone dihedrals
        vector<vector<double>> alpha;
        vector<vector<double>> beta;
        vector<vector<double>> gamma;
        vector<vector<double>> delta;
        vector<vector<double>> epsilon;
        vector<vector<double>> zeta;

        // Output: suite name (composed of a number-like character for the
        // 5' hemi-nucleotide and a letter-like character for the
        // 3' hemi-nucleotide) and suiteness score
        vector<string> suite_names;
        vector<string> suite_ddg;
        vector<double> suiteness;

        // Other internal variables
        size_t N_reference_ddg;
        vector<size_t> N_reference_suite;
        size_t N_continuous_group;
        vector<size_t> N_residue;
        size_t N_suite;
        double suiteness_cutoff;

    }; // RnaSuite class

}

#endif

