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

#include <RnaSuite.hpp>

using namespace std;

namespace loos {

    // |------------------------------------------------------------------------
    // | Constructors
    // |------------------------------------------------------------------------

    RnaSuite::RnaSuite(const AtomicGroup &group,
                       const double suiteness_cutoff_) {

        extractRnaBackboneAtoms(group);
        suiteness_cutoff = suiteness_cutoff_;

    }

    RnaSuite::RnaSuite(const AtomicGroup &group) {

        extractRnaBackboneAtoms(group);
        suiteness_cutoff = 0.01;

    }

    RnaSuite::RnaSuite() {
        suiteness_cutoff = 0.01;
    }

    // |------------------------------------------------------------------------
    // | Methods
    // |------------------------------------------------------------------------

    void RnaSuite::calculateBackboneDihedrals() {

        // Clear vector of vectors of doubles for each backbone dihedral
        alpha.clear();
        beta.clear();
        gamma.clear();
        delta.clear();
        epsilon.clear();
        zeta.clear();

        for (size_t i = 0; i < N_continuous_group; i++) {

            std::vector<double> continuous_alpha(N_residue[i]);
            std::vector<double> continuous_beta(N_residue[i]);
            std::vector<double> continuous_gamma(N_residue[i]);
            std::vector<double> continuous_delta(N_residue[i] + 1);
            std::vector<double> continuous_epsilon(N_residue[i]);
            std::vector<double> continuous_zeta(N_residue[i]);

            for (size_t j = 0; j < N_residue[i]; j++) {

                continuous_alpha[j] = Math::torsion(
                    alpha_atoms[i][j][0], alpha_atoms[i][j][1],
                    alpha_atoms[i][j][2], alpha_atoms[i][j][3]);
                continuous_beta[j] = Math::torsion(
                    beta_atoms[i][j][0], beta_atoms[i][j][1],
                    beta_atoms[i][j][2], beta_atoms[i][j][3]);
                continuous_gamma[j] = Math::torsion(
                    gamma_atoms[i][j][0], gamma_atoms[i][j][1],
                    gamma_atoms[i][j][2], gamma_atoms[i][j][3]);
                continuous_delta[j] = Math::torsion(
                    delta_atoms[i][j][0], delta_atoms[i][j][1],
                    delta_atoms[i][j][2], delta_atoms[i][j][3]);
                continuous_epsilon[j] = Math::torsion(
                    epsilon_atoms[i][j][0], epsilon_atoms[i][j][1],
                    epsilon_atoms[i][j][2], epsilon_atoms[i][j][3]);
                continuous_zeta[j] = Math::torsion(
                    zeta_atoms[i][j][0], zeta_atoms[i][j][1],
                    zeta_atoms[i][j][2], zeta_atoms[i][j][3]);

            }

            continuous_delta[N_residue[i]] = Math::torsion(
                delta_atoms[i][N_residue[i]][0],
                delta_atoms[i][N_residue[i]][1],
                delta_atoms[i][N_residue[i]][2],
                delta_atoms[i][N_residue[i]][3]);

            alpha.push_back(continuous_alpha);
            beta.push_back(continuous_beta);
            gamma.push_back(continuous_gamma);
            delta.push_back(continuous_delta);
            epsilon.push_back(continuous_epsilon);
            zeta.push_back(continuous_zeta);

        }

    } // calculateBackboneDihedrals()

    void RnaSuite::checkContinuousGroupSize(
        const std::vector<std::vector<AtomicGroup>> &group_vector,
        const size_t target_size, const string dihedral_name) const {

        if (group_vector.size() != target_size) {

            cout << boost::format("Error: different number of continuous "
                "groups for alpha (%d) and %s (%d)\n") % target_size
                % dihedral_name % group_vector.size();
            throw(LOOSError());

        }

    } // checkContinuousGroupSize()

    void RnaSuite::checkResidueSize(
        const std::vector<AtomicGroup> &residue_vector,
        const size_t target_size, const string dihedral_name,
        const size_t group_index) const {

        if (residue_vector.size() != target_size) {

            cout << boost::format("Error: different number of residues in "
                "continuous group %d for alpha (%d) and %s (%d)\n")
                % group_index % target_size % dihedral_name
                % residue_vector.size();

        }

    } // checkResidueSize()

    void RnaSuite::extractRnaBackboneAtoms(const AtomicGroup &group) {

        std::vector<AtomicGroup> continuous_alpha_atoms;
        std::vector<AtomicGroup> continuous_beta_atoms;
        std::vector<AtomicGroup> continuous_gamma_atoms;
        std::vector<AtomicGroup> continuous_delta_atoms;
        std::vector<AtomicGroup> continuous_epsilon_atoms;
        std::vector<AtomicGroup> continuous_zeta_atoms;
        AtomicGroup dihedral_atoms;
        AtomicGroup residue_p;
        AtomicGroup residue_o5p;
        AtomicGroup residue_c5p;
        AtomicGroup residue_c4p;
        AtomicGroup residue_c3p;
        AtomicGroup residue_o3p;
        AtomicGroup prev_residue_c4p;
        AtomicGroup prev_residue_c3p;
        AtomicGroup prev_residue_o3p;
        int current_resid = -2;

        // True if this is the initial residue in a continuous group
        bool first_res = true;

        // Clear vector of vectors of AtomicGroups for each backbone dihedral
        alpha_atoms.clear();
        beta_atoms.clear();
        gamma_atoms.clear();
        delta_atoms.clear();
        epsilon_atoms.clear();
        zeta_atoms.clear();

        // Extract all RNA backbone atoms (P, O5', C5', C4', C3', and O3') into
        // one AtomicGroup. Use raw string literal R"()" to avoid escaping "
        AtomicGroup backbone = selectAtoms(group,
                                           R"(name =~ "^(P|C[345]'|O[35]')$")");

        // Split by resid and loop over residues
        for (AtomicGroup residue : backbone.splitByResidue()) {

            // Select RNA backbone atoms from residue
            residue_p = selectAtoms(residue, R"(name == "P")");
            residue_o5p = selectAtoms(residue, R"(name == "O5'")");
            residue_c5p = selectAtoms(residue, R"(name == "C5'")");
            residue_c4p = selectAtoms(residue, R"(name == "C4'")");
            residue_c3p = selectAtoms(residue, R"(name == "C3'")");
            residue_o3p = selectAtoms(residue, R"(name == "O3'")");

            // If any atom besides P is missing, skip this residue and start a
            // new continuous group
            if (residue_o5p.size() != 1 || residue_c5p.size() != 1 ||
                residue_c4p.size() != 1 || residue_c3p.size() != 1 ||
                residue_o3p.size() != 1) {

                first_res = true;
                continue;

            }

            // If the resid is not sequential, this is not a continuous group
            if (residue_p.size() != 1
                || residue_p[0]->resid() != current_resid + 1) first_res = true;

            if (first_res) {

                first_res = false;

                // Record any previous continuous group
                if (continuous_alpha_atoms.size() != 0) {

                    alpha_atoms.push_back(continuous_alpha_atoms);
                    beta_atoms.push_back(continuous_beta_atoms);
                    gamma_atoms.push_back(continuous_gamma_atoms);
                    delta_atoms.push_back(continuous_delta_atoms);
                    epsilon_atoms.push_back(continuous_epsilon_atoms);
                    zeta_atoms.push_back(continuous_zeta_atoms);

                }

                // Clear vectors of AtomicGroups for this continuous groups
                continuous_alpha_atoms.clear();
                continuous_beta_atoms.clear();
                continuous_gamma_atoms.clear();
                continuous_delta_atoms.clear();
                continuous_epsilon_atoms.clear();
                continuous_zeta_atoms.clear();

                // Record delta for this initial residue
                dihedral_atoms = residue_c5p;
                dihedral_atoms.append(residue_c4p);
                dihedral_atoms.append(residue_c3p);
                dihedral_atoms.append(residue_o3p);
                continuous_delta_atoms.push_back(dihedral_atoms);

            } else {

                // Record backbone dihedrals for the remainder of the suite,
                // i.e. epsilon and zeta of the previous residue and alpha,
                // beta, gamma, and delta of the current residue
                dihedral_atoms = prev_residue_c4p;
                dihedral_atoms.append(prev_residue_c3p);
                dihedral_atoms.append(prev_residue_o3p);
                dihedral_atoms.append(residue_p);
                continuous_epsilon_atoms.push_back(dihedral_atoms);
                dihedral_atoms = prev_residue_c3p;
                dihedral_atoms.append(prev_residue_o3p);
                dihedral_atoms.append(residue_p);
                dihedral_atoms.append(residue_o5p);
                continuous_zeta_atoms.push_back(dihedral_atoms);
                dihedral_atoms = prev_residue_o3p;
                dihedral_atoms.append(residue_p);
                dihedral_atoms.append(residue_o5p);
                dihedral_atoms.append(residue_c5p);
                continuous_alpha_atoms.push_back(dihedral_atoms);
                dihedral_atoms = residue_p;
                dihedral_atoms.append(residue_o5p);
                dihedral_atoms.append(residue_c5p);
                dihedral_atoms.append(residue_c4p);
                continuous_beta_atoms.push_back(dihedral_atoms);
                dihedral_atoms = residue_o5p;
                dihedral_atoms.append(residue_c5p);
                dihedral_atoms.append(residue_c4p);
                dihedral_atoms.append(residue_c3p);
                continuous_gamma_atoms.push_back(dihedral_atoms);
                dihedral_atoms = residue_c5p;
                dihedral_atoms.append(residue_c4p);
                dihedral_atoms.append(residue_c3p);
                dihedral_atoms.append(residue_o3p);
                continuous_delta_atoms.push_back(dihedral_atoms);

            }

            // Save C4', C3', and O3' for dihedrals in the next residue
            prev_residue_c4p = residue_c4p;
            prev_residue_c3p = residue_c3p;
            prev_residue_o3p = residue_o3p;

            // Update resid
            current_resid = residue_o5p[0]->resid();

        } // loop over residues

        // Record any previous continuous group
        if (continuous_alpha_atoms.size() != 0) {

            alpha_atoms.push_back(continuous_alpha_atoms);
            beta_atoms.push_back(continuous_beta_atoms);
            gamma_atoms.push_back(continuous_gamma_atoms);
            delta_atoms.push_back(continuous_delta_atoms);
            epsilon_atoms.push_back(continuous_epsilon_atoms);
            zeta_atoms.push_back(continuous_zeta_atoms);

        }

        // Get number of continuous groups and check that all dihedral groups
        // have same size
        N_continuous_group = alpha_atoms.size();
        checkContinuousGroupSize(beta_atoms, N_continuous_group, "beta");
        checkContinuousGroupSize(gamma_atoms, N_continuous_group, "gamma");
        checkContinuousGroupSize(delta_atoms, N_continuous_group, "delta");
        checkContinuousGroupSize(epsilon_atoms, N_continuous_group, "epsilon");
        checkContinuousGroupSize(zeta_atoms, N_continuous_group, "zeta");

        // Get number of residues in each continuous group and check that these
        // are consistent across backbone dihedrals. Delta should have one
        // additional residue per continuous group.
        size_t residue_size;

        for (size_t i = 0; i < N_continuous_group; i++) {

            residue_size = alpha_atoms[i].size();
            checkResidueSize(beta_atoms[i], residue_size, "beta", i + 1);
            checkResidueSize(gamma_atoms[i], residue_size, "gamma", i + 1);
            checkResidueSize(delta_atoms[i], residue_size + 1, "delta", i + 1);
            checkResidueSize(epsilon_atoms[i], residue_size, "epsilon", i + 1);
            checkResidueSize(zeta_atoms[i], residue_size, "zeta", i + 1);
            N_residue.push_back(residue_size);

        }

    } // extractRnaBackboneAtoms()

    double RnaSuite::getSuitenessCutoff() const {
        return suiteness_cutoff;
    } // getSuitenessCutoff()

    void RnaSuite::printBackboneAtoms() const {

        size_t i_plus;
        size_t j_plus;

        cout << boost::format("Number of continuous groups: %d\n")
            % N_continuous_group;

        if (N_continuous_group == 0) return;

        for (size_t i = 0; i < N_continuous_group; i++) {

            i_plus = i + 1;
            cout << boost::format("Continuous group %d has %d residues\n")
                % i_plus % N_residue[i];

            for (size_t j = 0; j < N_residue[i]; j++) {

                j_plus = j + 1;
                cout << boost::format("Delta %d %d\n") % i_plus % j_plus;
                cout << delta_atoms[i][j] << endl;
                cout << boost::format("Epsilon %d %d\n") % i_plus % j_plus;
                cout << epsilon_atoms[i][j] << endl;
                cout << boost::format("Zeta %d %d\n") % i_plus % j_plus;
                cout << zeta_atoms[i][j] << endl;
                cout << boost::format("Alpha %d %d\n") % i_plus % j_plus;
                cout << alpha_atoms[i][j] << endl;
                cout << boost::format("Beta %d %d\n") % i_plus % j_plus;
                cout << beta_atoms[i][j] << endl;
                cout << boost::format("Gamma %d %d\n") % i_plus % j_plus;
                cout << gamma_atoms[i][j] << endl;

            }

            cout << boost::format("Delta %d %d\n") % i_plus
                % (N_residue[i] + 1);
            cout << delta_atoms[i][N_residue[i]] << endl;

        }

    } // printBackboneAtoms()

    void RnaSuite::printBackboneDihedrals() const {

        if (alpha.empty()) return;

        for (size_t i = 0; i < N_continuous_group; i++) {

            for (size_t j = 0; j < N_residue[i]; j++) {

                cout << boost::format("%4d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n")
                    % gamma_atoms[i][j][0]->resid() % delta[i][j]
                    % epsilon[i][j] % zeta[i][j] % alpha[i][j] % beta[i][j]
                    % gamma[i][j];

            }

            cout << boost::format("%4d %8.3f\n")
                % delta_atoms[i][N_residue[i]][0]->resid()
                % delta[i][N_residue[i]];

        }

    } // printBackboneDihedrals()

    void RnaSuite::setSuitenessCutoff(const double suiteness_cutoff_) {
        suiteness_cutoff = suiteness_cutoff_;
    } // setSuitenessCutoff()

}
