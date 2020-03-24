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

    // Constructors

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

    // Methods

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

    } // extractRnaBackboneAtoms()

    double RnaSuite::getSuitenessCutoff() const {
        return suiteness_cutoff;
    } // getSuitenessCutoff()

    void RnaSuite::printBackboneAtoms() const {

        uint continuous_counter;
        uint residue_counter;

        cout << boost::format("Sizes %d %d %d %d %d %d\n") % alpha_atoms.size()
            % beta_atoms.size() % gamma_atoms.size() % delta_atoms.size()
            % epsilon_atoms.size() % zeta_atoms.size();

        continuous_counter = 0;
        for (std::vector<AtomicGroup> continuous_atoms : alpha_atoms) {
            continuous_counter++;
            cout << boost::format("Alpha %d Size %d\n") % continuous_counter
                % continuous_atoms.size();
            residue_counter = 0;
            for (AtomicGroup residue_atoms : continuous_atoms) {
                residue_counter++;
                cout << boost::format("Alpha %d %d\n") % continuous_counter
                    % residue_counter;
                cout << residue_atoms << endl;
            }
        }

        continuous_counter = 0;
        for (std::vector<AtomicGroup> continuous_atoms : beta_atoms) {
            continuous_counter++;
            cout << boost::format("Beta %d Size %d\n") % continuous_counter
                % continuous_atoms.size();
            residue_counter = 0;
            for (AtomicGroup residue_atoms : continuous_atoms) {
                residue_counter++;
                cout << boost::format("Beta %d %d\n") % continuous_counter
                    % residue_counter;
                cout << residue_atoms << endl;
            }
        }

        continuous_counter = 0;
        for (std::vector<AtomicGroup> continuous_atoms : gamma_atoms) {
            continuous_counter++;
            cout << boost::format("Gamma %d Size %d\n") % continuous_counter
                % continuous_atoms.size();
            residue_counter = 0;
            for (AtomicGroup residue_atoms : continuous_atoms) {
                residue_counter++;
                cout << boost::format("Gamma %d %d\n") % continuous_counter
                    % residue_counter;
                cout << residue_atoms << endl;
            }
        }

        continuous_counter = 0;
        for (std::vector<AtomicGroup> continuous_atoms : delta_atoms) {
            continuous_counter++;
            cout << boost::format("Delta %d Size %d\n") % continuous_counter
                % continuous_atoms.size();
            residue_counter = 0;
            for (AtomicGroup residue_atoms : continuous_atoms) {
                residue_counter++;
                cout << boost::format("Delta %d %d\n") % continuous_counter
                    % residue_counter;
                cout << residue_atoms << endl;
            }
        }

        continuous_counter = 0;
        for (std::vector<AtomicGroup> continuous_atoms : epsilon_atoms) {
            continuous_counter++;
            cout << boost::format("Epsilon %d Size %d\n") % continuous_counter
                % continuous_atoms.size();
            residue_counter = 0;
            for (AtomicGroup residue_atoms : continuous_atoms) {
                residue_counter++;
                cout << boost::format("Epsilon %d %d\n") % continuous_counter
                    % residue_counter;
                cout << residue_atoms << endl;
            }
        }

        continuous_counter = 0;
        for (std::vector<AtomicGroup> continuous_atoms : zeta_atoms) {
            continuous_counter++;
            cout << boost::format("Zeta %d Size %d\n") % continuous_counter
                % continuous_atoms.size();
            residue_counter = 0;
            for (AtomicGroup residue_atoms : continuous_atoms) {
                residue_counter++;
                cout << boost::format("Zeta %d %d\n") % continuous_counter
                    % residue_counter;
                cout << residue_atoms << endl;
            }
        }

    } // printBackboneAtoms()

    void RnaSuite::setSuitenessCutoff(const double suiteness_cutoff_) {
        suiteness_cutoff = suiteness_cutoff_;
    } // setSuitenessCutoff()

}
