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
     *  This class acts on an AtomicGroup and assigns backbone suites to any RNA
     *  residues present. It also calculates the "suiteness" score that
     *  describes how well the residue fits into its assigned suite. The
     *  constructor requires that the user specifies a path to a file defining
     *  reference suites. The suites from Richardson et al. (2008) RNA 14,
     *  465-481 are included in $LOOS/share/suitename_definitions.dat
     */
    class RnaSuite {

    public:

        RnaSuite(const AtomicGroup &group, const string suite_defintion,
            const double suiteness_cutoff_);
        
        RnaSuite(const AtomicGroup &group, const string suite_definition);
        
        RnaSuite();
        
        //! Method to assign residues to backbone suites from Richardson et al.
        /**
         *  This method assigns residues to one of the reference suites defined
         *  in the constructor. The suite of a residue is defined from delta of
         *  the previous residue to delta of the current residue.
         */
        void assignSuitenameSuites();

        //! Method to calculate backbone dihedrals for each RNA residue
        /**
         *  This method calculates the six RNA backbone dihedrals (i.e. alpha,
         *  beta, gamma, delta, epsilon, and zeta) for each residue.
         */
        void calculateBackboneDihedrals();

        //! Method to define suites used for assignment
        /**
         *  This method defines reference suites. The argument must be a path to
         *  a file containing records consisting of fields with a width of eight
         *  characters. An example file for the suites defined in
         *  Richardson et al. (2008) RNA 14, 465-481 is included in
         *  $LOOS/share/suitename_definitions.dat. Records can be:
         *
         *  suite name ddg delta(i-1) epsilon zeta alpha beta gamma delta(i)
         *      Define a reference suite with name given in field 2, ddg_index
         *      given in field 3, and dihedrals of the cluster center given in
         *      fields 4 through 10.
         *
         *  width delta(i-1) epsilon zeta alpha beta gamma delta
         *      Define default widths for scaled hyperellipsoid distances.
         *
         *  domsat sat_name dom_name dihedral_index sat_width dom_width
         *      Define dominant-satellite pair with name of satellite suite in
         *      field 2, name of dominant suite in field 3, index of dihedral
         *      dimension with altered width in field 4, width of that dimension
         *      for satellite suite in field 5, and width of that dimension for
         *      dominant suite in field 6. Additional dimensions and width can
         *      be specified in fields 7 through 9, fields 10 through 12, etc.
         *
         *  dihedral min max
         *      Define allowed ranges for a dihedral. "dihedral" can be one of
         *      "delta", "epsilon", "zeta", "alpha", "beta", or "gamma". The
         *      minimum value is given in field 2 and maximum value in field 3.
         */
        void defineSuites(const string& suite_definition);

        //! Method to extract RNA backbone atoms from an AtomicGroup
        /**
         *  This method selects RNA backbone atoms (i.e. P, O5', C5', C4', C3',
         *  and O3') and splits them into AtomicGroups by residue id.
         */
        void extractRnaBackboneAtoms(const AtomicGroup &group);

        //! Method to return the current indices into delta delta gamma clusters
        vector<string> getSuiteDDGs() const;

        //! Method to return the current backbone dihedrals
        vector<vector<double> > getSuiteDihedrals() const;

        //! Method to return the current assigned suite names
        vector<string> getSuiteNames() const;

        //! Method to return the suite residue indices
        vector<int> getSuiteResids() const;

        //! Method to return the suite residue names
        vector<string> getSuiteResnames() const;

        //! Method to return the cutoff for the suiteness score of non-outliers
        double getSuitenessCutoff() const;

        //! Method to return the current suiteness scores
        vector<double> getSuitenessScores() const;

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

        //! Method to assign residues to a delta(i-1), delta, gamma index
        size_t assignDDGIndex(double dihedral, vector<double> &min,
            vector<double> &max, uint increment, uint &ddg_index);

        //! Calculate a dihedral in deg from 4 atoms in the range [0, 360]
        double calculateDihedral(const AtomicGroup &group);

        //! Method to check the size of a vector of continuous groups
        void checkContinuousGroupSize(
            const vector<vector<AtomicGroup> > &group_vector,
            const size_t target_size, const string dihedral_name) const;

        //! Method to check the size of a vector of residues
        void checkResidueSize(const vector<AtomicGroup> &residue_vector,
            const size_t target_size, const string dihedral_name,
            const size_t group_index) const;

        //! Method to test whether a point is in between two reference points
        bool isBetweenDomSatPair(const vector<double> &dihedrals,
            const vector<double> &dominant, const vector<double> &satellite);

        //! Calculate a scaled hyperellipsoid distance between two points
        double hyperellipsoidDist(const vector<double> &dihedrals,
            const vector<double> &reference, const vector<double> &width,
            uint first_index, uint last_index);

        // Reference suites used for assignment
        vector<vector<vector<double> >> reference_dihedrals;
        vector<vector<string> > reference_names;
        vector<string> reference_ddgs;

        // Widths used to scale each dihedral dimension
        vector<double> dihedral_width;

        // Indices of dominant-satellite pairs
        vector<vector<size_t> > dominant_suites;

        // Index into dominant-satellite pair widths
        vector<vector<size_t> > dom_sat_pair_index;

        // Alternative widths used to scale dominant-satellite pairs
        vector<vector<double> > dominant_width;
        vector<vector<double> > satellite_width;

        // Boundaries for allowed regions of delta(i-1), delta, and gamma
        vector<double> delta_min;
        vector<double> delta_max;
        vector<double> gamma_min;
        vector<double> gamma_max;

        // Boundaries for allowed regions of epsilon, zeta, alpha, beta
        vector<double> ezab_min;
        vector<double> ezab_max;

        // Vector of continuous groups, composed of vectors of AtomicGroups
        // for each residue within a continuous group
        vector<vector<AtomicGroup> > alpha_atoms;
        vector<vector<AtomicGroup> > beta_atoms;
        vector<vector<AtomicGroup> > gamma_atoms;
        vector<vector<AtomicGroup> > delta_atoms;
        vector<vector<AtomicGroup> > epsilon_atoms;
        vector<vector<AtomicGroup> > zeta_atoms;

        // Suite residue ids, residue names, and dihedrals
        vector<int> suite_resids;
        vector<string> suite_resnames;
        vector<vector<double> > suite_dihedrals;

        // Assigned suite names, ddg indices, and suiteness scores
        vector<string> suite_names;
        vector<string> suite_ddgs;
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

