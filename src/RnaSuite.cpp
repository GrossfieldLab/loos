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

    RnaSuite::RnaSuite(const AtomicGroup &group, const string suite_definition,
        const double suiteness_cutoff_) {

        suiteness_cutoff = suiteness_cutoff_;
        defineSuites(suite_definition);
        extractRnaBackboneAtoms(group);

    }

    RnaSuite::RnaSuite(const AtomicGroup &group,
        const string suite_definition) {

        suiteness_cutoff = 0.01;
        defineSuites(suite_definition);
        extractRnaBackboneAtoms(group);

    }

    RnaSuite::RnaSuite(const AtomicGroup &group,
        const double suiteness_cutoff_) {

        suiteness_cutoff = suiteness_cutoff_;
        defineSuites("suitename");
        extractRnaBackboneAtoms(group);

    }

    RnaSuite::RnaSuite(const AtomicGroup &group) {

        suiteness_cutoff = 0.01;
        defineSuites("suitename");
        extractRnaBackboneAtoms(group);

    }

    RnaSuite::RnaSuite() {

        suiteness_cutoff = 0.01;
        defineSuites("suitename");

    }

    // |------------------------------------------------------------------------
    // | Methods
    // |------------------------------------------------------------------------

    size_t RnaSuite::assignDDGIndex(double dihedral, vector<double> &min,
        vector<double> &max, uint increment, uint &ddg_index) {

        size_t i = 0;
        while (i < min.size()) {

            if (dihedral >= min[i] && dihedral <= max[i]) {

                ddg_index += i * increment;
                return i;

            }

            ++i;

        }

        return i;

    } // assignDDGIndex()

    void RnaSuite::assignSuitenameSuites() {

        size_t N_delta = delta_min.size();
        size_t N_gamma = gamma_min.size();
        size_t N_dg = N_delta * N_gamma;
        vector<double> suite(7);

        // Index into delta(i-1), delta, gamma clusters
        uint ddg_index;

        // Scaled 4D hyperellipsoid distance in epsilon, zeta, alpha, beta
        double dist_ezab;

        // Closest scaled 4D hyperellipsoid distance to a cluster and index of
        // the associated cluster
        double min_dist_ezab;
        size_t min_index;

        // Closest scaled 4D hyperellipsoid distance to a dominant cluster
        double dom_min_dist_ezab;
        size_t dom_min_index;

        // Closest scaled 4D hyperellipsoid distance to a non-dominant cluster
        double sat_min_dist_ezab;
        size_t sat_min_index;

        // Index into vector of widths for pair of dominant-satellite clusters
        size_t dom_sat_index;

        // Number of clusters this dinucleotide could belong to
        uint candidates;

        // Scaled 7D hyperellipsoid distance
        double dist_7;

        // Index of the assigned suite
        size_t assigned_suite_index;

        // Goodness-of-fit for assigned suite
        double suiteness_score;

        if (suite_dihedrals.empty()) {

            cerr << "Warning: backbone dihedrals are empty" << endl;
            return;

        }

        // Initialize vectors of suite names and suiteness scores
        suite_names.clear();
        suite_ddgs.clear();
        suiteness.clear();
        suite_names.reserve(N_suite);
        suite_ddgs.reserve(N_suite);
        suiteness.reserve(N_suite);

        for (size_t i = 0; i < N_suite; ++i) {

            // Assign delta(j-1), delta, gamma index. These 3 dihedrals have
            // 12 clusters that are independent of the other 4 dihedrals.
            ddg_index = 0;

            // Filter on 5' delta. Values outside of this range are
            // indicative of incorrect stereochemistry in the ribose.
            if (assignDDGIndex(suite_dihedrals[i][0], delta_min, delta_max,
                N_dg, ddg_index) == N_delta) {

                suite_names.push_back("!d");
                suite_ddgs.push_back("!!!");
                suiteness.push_back(0.0);
                continue;

            }

            // Filter on 3' delta
            if (assignDDGIndex(suite_dihedrals[i][6], delta_min, delta_max,
                N_gamma, ddg_index) == N_delta) {

                suite_names.push_back("!d");
                suite_ddgs.push_back("!!!");
                suiteness.push_back(0.0);
                continue;

            }

            // Filter on gamma
            if (assignDDGIndex(suite_dihedrals[i][5], gamma_min, gamma_max, 1,
                ddg_index) == N_gamma) {

                suite_names.push_back("!g");
                suite_ddgs.push_back("!!!");
                suiteness.push_back(0.0);
                continue;

            }

            // Filter on epsilon. Values outside of this range are
            // indicative of a misfit sugar pucker.
            if (suite_dihedrals[i][1] < ezab_min[0]
                || suite_dihedrals[i][1] > ezab_max[0]) {

                suite_names.push_back("!e");
                suite_ddgs.push_back("!!!");
                suiteness.push_back(0.0);
                continue;

            }

            // Filter on zeta
            if (suite_dihedrals[i][2] < ezab_min[1]
                || suite_dihedrals[i][2] > ezab_max[1]) {

                suite_names.push_back("!z");
                suite_ddgs.push_back("!!!");
                suiteness.push_back(0.0);
                continue;

            }

            // Filter on alpha
            if (suite_dihedrals[i][3] < ezab_min[2]
                || suite_dihedrals[i][3] > ezab_max[2]) {

                suite_names.push_back("!a");
                suite_ddgs.push_back("!!!");
                suiteness.push_back(0.0);
                continue;

            }

            // Filter on beta
            if (suite_dihedrals[i][4] < ezab_min[3]
                || suite_dihedrals[i][4] > ezab_max[3]) {

                suite_names.push_back("!b");
                suite_ddgs.push_back("!!!");
                suiteness.push_back(0.0);
                continue;

            }

            // If there are no clusters associated with this ddg_index, then
            // this is an outlier
            if (N_reference_suite[ddg_index] == 0) {

                suite_names.push_back("!!");
                suite_ddgs.push_back(reference_ddgs[ddg_index]);
                suiteness.push_back(0.0);
                continue;

            }

            // Find closest cluster in epsilon, zeta, alpha, beta
            // Largest distance in 7D is 688.66^3, so 10^9 should be safe
            min_dist_ezab = 999999999.0;
            dom_min_dist_ezab = 999999999.0;
            sat_min_dist_ezab = 999999999.0;
            min_index = N_reference_suite[ddg_index];
            dom_min_index = N_reference_suite[ddg_index];
            sat_min_index = N_reference_suite[ddg_index];
            candidates = 0;

            for (size_t j = 0; j < N_reference_suite[ddg_index]; ++j) {

                // Get 4D scaled hyperellipsoid distance
                dist_ezab = hyperellipsoidDist(suite_dihedrals[i],
                    reference_dihedrals[ddg_index][j], dihedral_width, 1, 4);

                // Get closest cluster
                if (dist_ezab < min_dist_ezab) {

                    min_dist_ezab = dist_ezab;
                    min_index = j;

                }

                // Get closest non-dominant cluster
                if (dominant_suites[ddg_index][j] != j
                    && dist_ezab < sat_min_dist_ezab) {

                    sat_min_dist_ezab = dist_ezab;
                    sat_min_index = j;

                }

                // If 4D distance < 1, this reference suite is a candidate
                if (dist_ezab < 1) {

                    ++candidates;

                    // Is this candidate a dominant cluster?
                    if (dominant_suites[ddg_index][j] == j) {

                        dom_min_dist_ezab = dist_ezab;
                        dom_min_index = j;

                    }

                }

            } // loop over reference suites

            // Assign membership to a reference suite

            // If there are multiple candidates, and the two canidates are
            // a dominant-satellite pair, then reweight distances
            if (candidates > 1 && dom_min_index != N_reference_suite[ddg_index]
                && sat_min_index != N_reference_suite[ddg_index]
                && dominant_suites[ddg_index][sat_min_index] == dom_min_index) {

                // Is the DNMP in between the dominant and satellite suites?
                if (isBetweenDomSatPair(suite_dihedrals[i],
                        reference_dihedrals[ddg_index][dom_min_index],
                        reference_dihedrals[ddg_index][sat_min_index])) {

                    // Rescale distances from point to dominant and satellite
                    // suites by ratio of distances from suite centers to
                    // boundary plane and assign to closest of the two
                    dom_sat_index = dom_sat_pair_index[ddg_index][sat_min_index];
                    if (hyperellipsoidDist(suite_dihedrals[i],
                            reference_dihedrals[ddg_index][sat_min_index],
                            satellite_width[dom_sat_index], 1, 4)
                        <= hyperellipsoidDist(suite_dihedrals[i],
                            reference_dihedrals[ddg_index][dom_min_index],
                            dominant_width[dom_sat_index], 1, 4))

                        assigned_suite_index = sat_min_index;

                    else assigned_suite_index = dom_min_index;
                        

                }

                else {

                    // Assign to closer of dominant or satellite suite
                    if (sat_min_dist_ezab <= dom_min_dist_ezab)
                        assigned_suite_index = sat_min_index;
                    else assigned_suite_index = dom_min_index;

                }

            }

            // If there is zero or one candidate or multiple candidates but no
            // dominant-satellite pair, then assign to the closest suite
            else assigned_suite_index = min_index;

            // Make a final decision on whether this is an outlier using 7D
            // hyperellipsoid distance
            dist_7 = hyperellipsoidDist(suite_dihedrals[i],
                reference_dihedrals[ddg_index][assigned_suite_index],
                dihedral_width, 0, 6);

            if (dist_7 < 1) {

                suite_names.push_back(
                    reference_names[ddg_index][assigned_suite_index]);
                suite_ddgs.push_back(reference_ddgs[ddg_index]);
                suiteness_score = (1 + cos(M_PI * cbrt(dist_7))) / 2.0;
                if (suiteness_score < suiteness_cutoff)
                    suiteness_score = suiteness_cutoff;
                suiteness.push_back(suiteness_score);

            } else {

                suite_names.push_back("!!");
                suite_ddgs.push_back(reference_ddgs[ddg_index]);
                suiteness.push_back(0.0);

            }

        } // loop over suites

    } // assignSuitenameSuites()

    void RnaSuite::calculateBackboneDihedrals() {

        double prev_delta;
        vector<double> suite(7);

        // Clear vector of doubles for suite backbone dihedrals
        suite_dihedrals.clear();

        for (size_t i = 0; i < N_continuous_group; ++i) {

            prev_delta = calculateDihedral(delta_atoms[i][0]);

            for (size_t j = 0; j < N_residue[i]; ++j) {

                suite[0] = prev_delta;
                suite[1] = calculateDihedral(epsilon_atoms[i][j]);
                suite[2] = calculateDihedral(zeta_atoms[i][j]);
                suite[3] = calculateDihedral(alpha_atoms[i][j]);
                suite[4] = calculateDihedral(beta_atoms[i][j]);
                suite[5] = calculateDihedral(gamma_atoms[i][j]);
                prev_delta = calculateDihedral(delta_atoms[i][j + 1]);
                suite[6] = prev_delta;
                suite_dihedrals.push_back(suite);

            }

        }

    } // calculateBackboneDihedrals()

    double RnaSuite::calculateDihedral(const AtomicGroup &group) {

        double dihedral = Math::torsion(group[0], group[1], group[2], group[3]);
        if (dihedral < 0.0) dihedral += 360.0;
        return dihedral;

    } // calculateDihedral()

    void RnaSuite::checkContinuousGroupSize(
        const vector<vector<AtomicGroup> > &group_vector,
        const size_t target_size, const string dihedral_name) const {

        if (group_vector.size() != target_size) {

            cerr << boost::format("Error: different number of continuous "
                "groups for alpha (%d) and %s (%d)") % target_size
                % dihedral_name % group_vector.size() << endl;
            throw(LOOSError());

        }

    } // checkContinuousGroupSize()

    void RnaSuite::checkResidueSize(
        const vector<AtomicGroup> &residue_vector,
        const size_t target_size, const string dihedral_name,
        const size_t group_index) const {

        if (residue_vector.size() != target_size) {

            cerr << boost::format("Error: different number of residues in "
                "continuous group %d for alpha (%d) and %s (%d)") % group_index
                % target_size % dihedral_name % residue_vector.size() << endl;
            throw(LOOSError());

        }

    } // checkResidueSize()

    void RnaSuite::defineSuites(const string& suite_definition) {

        // Clear vectors for reference suites
        reference_dihedrals.clear();
        reference_names.clear();
        reference_ddgs.clear();
        dihedral_width.clear();
        dominant_suites.clear();
        dom_sat_pair_index.clear();
        dominant_width.clear();
        satellite_width.clear();
        delta_min.clear();
        delta_max.clear();
        gamma_min.clear();
        gamma_max.clear();
        ezab_min = vector<double>(4);
        ezab_max = vector<double>(4);
        N_reference_ddg = 0;
        N_reference_suite.clear();

        if (suite_definition == "suitename")
            defineSuitesFromFile("suitename_definitions.dat");

        else defineSuitesFromFile(suite_definition);

    } // defineSuites()

    void RnaSuite::defineSuitesFromFile(const string& filename) {

        size_t ddg_index;
        size_t dom_index;
        size_t sat_index;
        size_t position;
        string field;
        string line;
        string record;
        vector<double> dihedrals(7);

        // Store dominant-satellite pairs
        vector<size_t> domsat_ddg;
        vector<size_t> domsat_dom;
        vector<size_t> domsat_sat;
        vector<vector<size_t> > domsat_dihedral;
        vector<vector<double> > domsat_dom_width;
        vector<vector<double> > domsat_sat_width;

        // Read file contents
        ifstream ifs(filename.c_str());
        if (!ifs) throw(FileOpenError(filename));

        while (getline(ifs, line)) {

            record = parseStringAs<std::string>(line, 0, 8);

            if (record.empty() || record[0] == '#') continue;

            else if (record == "suite") {

                // Define a reference suite

                // Get delta delta gamma cluster
                field = parseStringAs<std::string>(
                    line, 16, min((size_t) 8, line.size() - 16));
                ddg_index = N_reference_ddg;
                for (size_t i = 0; i < N_reference_ddg; ++i)
                    if (field == reference_ddgs[i]) {
                        ddg_index = i;
                        break;
                    }

                // This is a new DDG cluster
                if (ddg_index == N_reference_ddg) {
                    reference_ddgs.push_back(field);
                    reference_dihedrals.push_back(vector<vector<double> >());
                    reference_names.push_back(vector<string>());
                    ++N_reference_ddg;
                    N_reference_suite.push_back(0);
                }

                // Get suite name
                field = parseStringAs<std::string>(line, 8, 8);
                if (field.empty()) continue;
                reference_names[ddg_index].push_back(field);

                // Get reference suite dihedrals
                dihedrals[0] = parseStringAs<double>(line, 24, 8);
                dihedrals[1] = parseStringAs<double>(line, 32, 8);
                dihedrals[2] = parseStringAs<double>(line, 40, 8);
                dihedrals[3] = parseStringAs<double>(line, 48, 8);
                dihedrals[4] = parseStringAs<double>(line, 56, 8);
                dihedrals[5] = parseStringAs<double>(line, 64, 8);
                dihedrals[6] = parseStringAs<double>(line, 72, 8);
                reference_dihedrals[ddg_index].push_back(dihedrals);

                ++N_reference_suite[ddg_index];

            } else if (record == "width") {

                // Get default widths for hyperellipsoid distance
                dihedral_width.push_back(parseStringAs<double>(line, 8, 8));
                dihedral_width.push_back(parseStringAs<double>(line, 16, 8));
                dihedral_width.push_back(parseStringAs<double>(line, 24, 8));
                dihedral_width.push_back(parseStringAs<double>(line, 32, 8));
                dihedral_width.push_back(parseStringAs<double>(line, 40, 8));
                dihedral_width.push_back(parseStringAs<double>(line, 48, 8));
                dihedral_width.push_back(parseStringAs<double>(line, 56, 8));

            } else if (record == "domsat") {

                // Define a dominant-satellite pair

                // Get index of dominant suite
                field = parseStringAs<std::string>(line, 16, 8);
                ddg_index = N_reference_ddg;
                for (size_t i = 0; i < N_reference_ddg; ++i) {
                    for (size_t j = 0; j < N_reference_suite[i]; ++j)
                        if (field == reference_names[i][j]) {
                            ddg_index = i;
                            dom_index = j;
                            break;
                        }
                    if (ddg_index != N_reference_ddg) break;
                }

                if (ddg_index == N_reference_ddg) {
                    cerr << boost::format(
                        "Warning: dominant suite %s was not defined in file %s")
                        % field % filename << endl;
                    continue;
                }

                // Get index of satellite suite
                field = parseStringAs<std::string>(line, 8, 8);
                sat_index = N_reference_suite[ddg_index];
                for (size_t j = 0; j < N_reference_suite[ddg_index]; ++j)
                    if (field == reference_names[ddg_index][j]) {
                        sat_index = j;
                        break;
                    }

                if (sat_index == N_reference_suite[ddg_index]) {
                    cerr << boost::format(
                        "Warning: satellite suite %s was not defined in file %s")
                        % field % filename << endl;
                    continue;
                }

                domsat_ddg.push_back(ddg_index);
                domsat_dom.push_back(dom_index);
                domsat_sat.push_back(sat_index);

                // Loop over dihedrals with alternate widths
                vector<size_t> dihedral_indices;
                vector<double> dom_width;
                vector<double> sat_width;
                position = 24;
                while (position < line.size()) {
                    dihedral_indices.push_back(
                        parseStringAs<size_t>(line, position, 8));
                    sat_width.push_back(
                        parseStringAs<double>(line, position + 8, 8));
                    dom_width.push_back(
                        parseStringAs<double>(line, position + 16, 8));
                    position += 24;
                }
                domsat_dihedral.push_back(dihedral_indices);
                domsat_dom_width.push_back(dom_width);
                domsat_sat_width.push_back(sat_width);

            } else if (record == "delta") {

                delta_min.push_back(parseStringAs<double>(line, 8, 8));
                delta_max.push_back(parseStringAs<double>(line, 16, 8));

            } else if (record == "epsilon") {

                ezab_min[0] = parseStringAs<double>(line, 8, 8);
                ezab_max[0] = parseStringAs<double>(line, 16, 8);

            } else if (record == "zeta") {

                ezab_min[1] = parseStringAs<double>(line, 8, 8);
                ezab_max[1] = parseStringAs<double>(line, 16, 8);

            } else if (record == "alpha") {

                ezab_min[2] = parseStringAs<double>(line, 8, 8);
                ezab_max[2] = parseStringAs<double>(line, 16, 8);

            } else if (record == "beta") {

                ezab_min[3] = parseStringAs<double>(line, 8, 8);
                ezab_max[3] = parseStringAs<double>(line, 16, 8);

            } else if (record == "gamma") {

                gamma_min.push_back(parseStringAs<double>(line, 8, 8));
                gamma_max.push_back(parseStringAs<double>(line, 16, 8));

            } else cerr << boost::format(
                "Warning: Unrecognized record %s in suite definition from %s")
                % record % filename << endl;

        } // Loop over lines in file

        // Construct vectors for dominant-satellite pairs
        for (size_t i = 0; i < N_reference_ddg; ++i) {
            dominant_suites.push_back(
                vector<size_t>(N_reference_suite[i], N_reference_suite[i]));
            dom_sat_pair_index.push_back(
                vector<size_t>(N_reference_suite[i], domsat_dihedral.size()));
        }

        for (size_t i = 0; i < domsat_dihedral.size(); ++i) {
            dominant_suites[domsat_ddg[i]][domsat_dom[i]] = domsat_dom[i];
            dominant_suites[domsat_ddg[i]][domsat_sat[i]] = domsat_dom[i];
            dom_sat_pair_index[domsat_ddg[i]][domsat_sat[i]] = i;
            dominant_width.push_back(dihedral_width);
            satellite_width.push_back(dihedral_width);
            for (size_t j = 0; j < domsat_dihedral[i].size(); ++j) {
                dominant_width[i][domsat_dihedral[i][j]] = domsat_dom_width[i][j];
                satellite_width[i][domsat_dihedral[i][j]] = domsat_sat_width[i][j];
            }
        }

    } // defineSuitesFromFile()

    void RnaSuite::extractRnaBackboneAtoms(const AtomicGroup &group) {

        vector<AtomicGroup> continuous_alpha_atoms;
        vector<AtomicGroup> continuous_beta_atoms;
        vector<AtomicGroup> continuous_gamma_atoms;
        vector<AtomicGroup> continuous_delta_atoms;
        vector<AtomicGroup> continuous_epsilon_atoms;
        vector<AtomicGroup> continuous_zeta_atoms;
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
        size_t residue_size;

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
            "(name =~ \"^(P|C[345]'|O[35]')$\")");

        // Split by resid and loop over residues
        vector<AtomicGroup> backbone_residues = backbone.splitByResidue();
        for (size_t i = 0; i < backbone_residues.size(); ++i) {

            // Select RNA backbone atoms from residue
            residue_p = selectAtoms(backbone_residues[i], "(name == \"P\")");
            residue_o5p = selectAtoms(backbone_residues[i], "(name == \"O5'\")");
            residue_c5p = selectAtoms(backbone_residues[i], "(name == \"C5'\")");
            residue_c4p = selectAtoms(backbone_residues[i], "(name == \"C4'\")");
            residue_c3p = selectAtoms(backbone_residues[i], "(name == \"C3'\")");
            residue_o3p = selectAtoms(backbone_residues[i], "(name == \"O3'\")");

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
        N_residue.clear();
        suite_resids.clear();
        suite_resnames.clear();

        for (size_t i = 0; i < N_continuous_group; ++i) {

            residue_size = alpha_atoms[i].size();
            checkResidueSize(beta_atoms[i], residue_size, "beta", i + 1);
            checkResidueSize(gamma_atoms[i], residue_size, "gamma", i + 1);
            checkResidueSize(delta_atoms[i], residue_size + 1, "delta", i + 1);
            checkResidueSize(epsilon_atoms[i], residue_size, "epsilon", i + 1);
            checkResidueSize(zeta_atoms[i], residue_size, "zeta", i + 1);
            N_residue.push_back(residue_size);

            for (size_t j = 0; j < residue_size; ++j) {

                suite_resids.push_back(gamma_atoms[i][j][0]->resid());
                suite_resnames.push_back(gamma_atoms[i][j][0]->resname());

            }

        }

        N_suite = suite_resids.size();

    } // extractRnaBackboneAtoms()

    vector<string> RnaSuite::getSuiteDDGs() const {
        return suite_ddgs;
    } // getSuiteDDGs()

    vector<vector<double> > RnaSuite::getSuiteDihedrals() const {
        return suite_dihedrals;
    } // getSuiteDihedrals()

    vector<string> RnaSuite::getSuiteNames() const {
        return suite_names;
    } // getSuiteNames()

    vector<int> RnaSuite::getSuiteResids() const {
        return suite_resids;
    } // getSuiteResids()

    vector<string> RnaSuite::getSuiteResnames() const {
        return suite_resnames;
    } // getSuiteResnames()

    double RnaSuite::getSuitenessCutoff() const {
        return suiteness_cutoff;
    } // getSuitenessCutoff()

    vector<double> RnaSuite::getSuitenessScores() const {
        return suiteness;
    } // getSuitenessScores()

    double RnaSuite::hyperellipsoidDist(const vector<double> &dihedrals,
        const vector<double> &reference, const vector<double> &width,
        uint first_index, uint last_index) {

        double unscaled_diff;
        double sum_scaled_powers = 0.0;

        for (uint i = first_index; i <= last_index; ++i) {

            unscaled_diff = abs(dihedrals[i] - reference[i]);
            // suitename program does not wrap unscaled coordinates
            // if (unscaled_diff > 180.0) unscaled_diff = 360.0 - unscaled_diff;
            sum_scaled_powers += pow(unscaled_diff / width[i], 3.0);

        }

        return sum_scaled_powers;

    } // hyperellipsoidDist4()

    bool RnaSuite::isBetweenDomSatPair(const vector<double> &dihedrals,
        const vector<double> &dominant, const vector<double> &satellite) {

        double dom_to_sat;
        double dom_dot_product = 0;
        double sat_dot_product = 0;

        // If the point is in between the dominant and satellite reference
        // suites, then the dot product between the vectors (point - dominant)
        // and (satellite - dominant) and the dot product between the vectors
        // (point - satellite) and (dominant - satellite) should both be
        // positive, i.e. the cosine of the angles is positive.
        for (uint i = 1; i <= 4; ++i) {

            dom_to_sat = satellite[i] - dominant[i];
            dom_dot_product += (dihedrals[i] - dominant[i]) * dom_to_sat;
            // sat_dot_product += (dihedrals[i] - satellite[i]) * sat_to_dom
            // sat_to_dom = -dom_to_sat
            sat_dot_product += (satellite[i] - dihedrals[i]) * dom_to_sat;

        }

        return dom_dot_product > 0 && sat_dot_product > 0;

    } // isBetweenDomSatPair()

    void RnaSuite::printBackboneAtoms() const {

        size_t i_plus;
        size_t j_plus;

        cout << "\n    ====  Printing backbone atoms  ====\n" << endl;

        if (alpha_atoms.empty()) {

            cerr << "Warning: backbone atoms are empty" << endl;
            return;

        }

        cout << boost::format("Number of continuous groups: %d\n")
            % N_continuous_group;

        for (size_t i = 0; i < N_continuous_group; ++i) {

            i_plus = i + 1;
            cout << boost::format("Continuous group %d has %d residues\n")
                % i_plus % N_residue[i];

            for (size_t j = 0; j < N_residue[i]; ++j) {

                j_plus = j + 1;
                cout << boost::format("Delta %d %d\n") % i_plus % j_plus;
                cout << delta_atoms[i][j] << "\n";
                cout << boost::format("Epsilon %d %d\n") % i_plus % j_plus;
                cout << epsilon_atoms[i][j] << "\n";
                cout << boost::format("Zeta %d %d\n") % i_plus % j_plus;
                cout << zeta_atoms[i][j] << "\n";
                cout << boost::format("Alpha %d %d\n") % i_plus % j_plus;
                cout << alpha_atoms[i][j] << "\n";
                cout << boost::format("Beta %d %d\n") % i_plus % j_plus;
                cout << beta_atoms[i][j] << "\n";
                cout << boost::format("Gamma %d %d\n") % i_plus % j_plus;
                cout << gamma_atoms[i][j] << "\n";

            }

            cout << boost::format("Delta %d %d\n") % i_plus
                % (N_residue[i] + 1);
            cout << delta_atoms[i][N_residue[i]] << "\n";

        }

    } // printBackboneAtoms()

    void RnaSuite::printBackboneDihedrals() const {

        cout << "\n    ====  Printing backbone dihedrals  ====\n" << endl;

        if (suite_dihedrals.empty()) {

            cerr << "Warning: backbone dihedrals are empty" << endl;
            return;

        }

        for (size_t i = 0; i < N_suite; ++i)
            cout << boost::format(
                "%5d %3s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n")
                % suite_resids[i] % suite_resnames[i] % suite_dihedrals[i][0]
                % suite_dihedrals[i][1] % suite_dihedrals[i][2]
                % suite_dihedrals[i][3] % suite_dihedrals[i][4]
                % suite_dihedrals[i][5] % suite_dihedrals[i][6];

    } // printBackboneDihedrals()

    void RnaSuite::printReferenceSuites() const {

        cout << "\n    ====  Printing reference suites  ====\n" << endl;

        if (reference_dihedrals.empty()) {

            cerr << "Warning: reference suites are empty" << endl;
            return;

        }

        for (size_t i = 0; i < N_reference_ddg; ++i)
            for (size_t j = 0; j < N_reference_suite[i]; ++j)
                cout << boost::format(
                    "%2s %3s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n")
                    % reference_names[i][j] % reference_ddgs[i]
                    % reference_dihedrals[i][j][0]
                    % reference_dihedrals[i][j][1]
                    % reference_dihedrals[i][j][2]
                    % reference_dihedrals[i][j][3]
                    % reference_dihedrals[i][j][4]
                    % reference_dihedrals[i][j][5]
                    % reference_dihedrals[i][j][6];

    } // printReferenceSuites()

    void RnaSuite::printSuites() const {

        cout << "\n    ====  Printing suites  ====\n" << endl;

        if (suite_names.empty()) {

            cerr << "Warning: suites are empty" << endl;
            return;

        }

        for (size_t i = 0; i < N_suite; ++i)
            cout << boost::format("%5d %3s %2s %3s %8.6f %7.3f %7.3f %7.3f "
                "%7.3f %7.3f %7.3f %7.3f\n") % suite_resids[i]
                % suite_resnames[i] % suite_names[i] % suite_ddgs[i]
                % suiteness[i] % suite_dihedrals[i][0] % suite_dihedrals[i][1]
                % suite_dihedrals[i][2] % suite_dihedrals[i][3]
                % suite_dihedrals[i][4] % suite_dihedrals[i][5]
                % suite_dihedrals[i][6];

    } // printSuites()

    void RnaSuite::setSuitenessCutoff(const double suiteness_cutoff_) {

        suiteness_cutoff = suiteness_cutoff_;

    } // setSuitenessCutoff()

}

