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
            if (suite_dihedrals[i][1] < filter_min[0]
                || suite_dihedrals[i][1] > filter_max[0]) {

                suite_names.push_back("!e");
                suite_ddgs.push_back("!!!");
                suiteness.push_back(0.0);
                continue;

            }

            // Filter on zeta
            if (suite_dihedrals[i][2] < filter_min[1]
                || suite_dihedrals[i][2] > filter_max[1]) {

                suite_names.push_back("!z");
                suite_ddgs.push_back("!!!");
                suiteness.push_back(0.0);
                continue;

            }

            // Filter on alpha
            if (suite_dihedrals[i][3] < filter_min[2]
                || suite_dihedrals[i][3] > filter_max[2]) {

                suite_names.push_back("!a");
                suite_ddgs.push_back("!!!");
                suiteness.push_back(0.0);
                continue;

            }

            // Filter on beta
            if (suite_dihedrals[i][4] < filter_min[3]
                || suite_dihedrals[i][4] > filter_max[3]) {

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
        const vector<vector<AtomicGroup>> &group_vector,
        const size_t target_size, const string dihedral_name) const {

        if (group_vector.size() != target_size) {

            cerr << boost::format("Error: different number of continuous "
                "groups for alpha (%d) and %s (%d)\n") % target_size
                % dihedral_name % group_vector.size();
            throw(LOOSError());

        }

    } // checkContinuousGroupSize()

    void RnaSuite::checkResidueSize(
        const vector<AtomicGroup> &residue_vector,
        const size_t target_size, const string dihedral_name,
        const size_t group_index) const {

        if (residue_vector.size() != target_size) {

            cerr << boost::format("Error: different number of residues in "
                "continuous group %d for alpha (%d) and %s (%d)\n")
                % group_index % target_size % dihedral_name
                % residue_vector.size();

        }

    } // checkResidueSize()

    void RnaSuite::defineSuites(const string suite_definition) {

        reference_dihedrals.clear();
        reference_names.clear();
        reference_ddgs.clear();

        if (suite_definition == "suitename"
            || suite_definition == "richardson") defineSuitesFromSuitename();

        else {

            cerr << boost::format("%s is not a recognized suite definition\n")
                % suite_definition;
            cerr << "Must be one of: suitename" << endl;
            throw(LOOSError());

        }

    } // defineSuites()

    void RnaSuite::defineSuitesFromFile(const string filename) {

        // TODO read suite definitions from file
        cerr << "Reading suite definitions from a file is not yet supported\n"
            "Go yell at Chapin" << endl;

    } // defineSuitesFromFile()

    void RnaSuite::defineSuitesFromSuitename() {

        size_t suite_size;

        // Means of dihedral angles
        reference_dihedrals = {
            { // ddg index 0: C3' C3' plus
                { 81.495, 212.250, 288.831, 294.967, 173.990,  53.550,  81.035},
                { 83.513, 218.120, 291.593, 292.247, 222.300,  58.067,  86.093},
                { 85.664, 245.014, 268.257, 303.879, 138.164,  61.950,  79.457},
                { 82.112, 190.682, 264.945, 295.967, 181.839,  51.455,  81.512},
                { 83.414, 217.400, 222.006, 302.856, 160.719,  49.097,  82.444},
                { 85.072, 216.324, 173.276, 289.320, 164.132,  45.876,  84.956},
                { 83.179, 210.347, 121.474, 288.568, 157.268,  49.347,  81.047},
                { 80.888, 218.636, 290.735, 167.447, 159.565,  51.326,  85.213},
                { 83.856, 238.750, 256.875,  69.562, 170.200,  52.800,  85.287},
                { 85.295, 244.085, 203.815,  65.880, 181.130,  54.680,  86.035},
                { 79.671, 202.471,  63.064,  68.164, 143.450,  49.664,  82.757},
                { 84.000, 195.000, 146.000, 170.000, 170.000,  52.000,  84.000}
            }, { // ddg index 1: C3' C3' trans
                { 80.514, 200.545, 280.510, 249.314,  82.662, 167.890,  85.507},
                { 80.223, 196.591, 291.299, 153.060, 194.379, 179.061,  83.648},
                { 81.395, 203.030, 294.445, 172.195, 138.540, 175.565,  84.470},
                { 87.417, 223.558,  80.175,  66.667, 109.150, 176.475,  83.833},
                { 86.055, 246.502, 100.392,  73.595, 213.752, 183.395,  85.483}
            }, { // ddg index 2: C3' C3' minus
            }, { // ddg index 3: C3' C2' plus
                { 84.215, 215.014, 288.672, 300.420, 177.476,  58.307, 144.841},
                { 82.731, 220.463, 288.665, 296.983, 221.654,  54.213, 143.771},
                { 84.700, 226.400, 168.336, 292.771, 177.629,  48.629, 147.950},
                { 83.358, 206.042, 277.567, 195.700, 161.600,  50.750, 145.258},
                { 82.614, 206.440,  52.524, 163.669, 148.421,  50.176, 147.590},
                { 84.285, 236.600, 220.400,  68.300, 200.122,  53.693, 145.730},
                { 84.457, 213.286,  69.086,  75.500, 156.671,  57.486, 147.686}
            }, { // ddg index 4: C3' C2' trans
                { 81.200, 199.243, 288.986, 180.286, 194.743, 178.200, 147.386},
                { 82.133, 204.933,  69.483,  63.417, 115.233, 176.283, 145.733}
            }, { // ddg index 5: C3' C2' minus
                { 83.977, 216.508, 287.192, 297.254, 225.154, 293.738, 150.677},
                { 84.606, 232.856, 248.125,  63.269, 181.975, 295.744, 149.744},
                { 83.000, 196.900,  65.350,  60.150, 138.425, 292.550, 154.275}
            }, { // ddg index 6: C2' C3' plus
                {145.399, 260.339, 288.756, 288.444, 192.733,  53.097,  84.067},
                {146.275, 259.783, 169.958, 298.450, 169.583,  50.908,  83.967},
                {149.286, 223.159, 139.421, 284.559, 158.107,  47.900,  84.424},
                {148.006, 191.944, 146.231, 289.288, 150.781,  42.419,  84.956},
                {148.028, 256.922, 165.194, 204.961, 165.194,  49.383,  82.983},
                {145.337, 262.869,  79.588, 203.863, 189.688,  58.000,  84.900},
                {148.992, 270.596, 240.892,  62.225, 176.271,  53.600,  87.262},
                {149.822, 249.956, 187.678,  80.433, 198.133,  61.000,  89.378},
                {146.922, 241.222,  88.894,  59.344, 160.683,  52.333,  83.417},
                {141.900, 258.383, 286.517, 178.267, 165.217,  48.350,  84.783}
            }, { // ddg index 7: C2' C3' trans
                {147.782, 260.712, 290.424, 296.200, 177.282, 175.594,  86.565},
                {143.722, 227.256, 203.789,  73.856, 216.733, 194.444,  80.911},
                {148.717, 274.683, 100.283,  80.600, 248.133, 181.817,  82.600},
                {150.311, 268.383,  84.972,  63.811, 191.483, 176.644,  85.600},
                {141.633, 244.100,  66.056,  71.667, 122.167, 182.200,  83.622}
            }, { // ddg index 8: C2' C3' minus
                {149.070, 249.780, 111.520, 278.370, 207.780, 287.820,  86.650}
            }, { // ddg index 9: C2' C2' plus
                {146.383, 259.402, 291.275, 291.982, 210.048,  54.412, 147.760},
                {145.256, 244.622, 162.822, 294.159, 171.630,  45.900, 145.804},
                {147.593, 248.421, 112.086, 274.943, 164.764,  56.843, 146.264},
                {150.077, 260.246, 213.785,  71.900, 207.638,  56.715, 148.131},
                {146.415, 257.831,  89.597,  67.923, 173.051,  55.513, 147.623},
                {142.900, 236.550, 268.800, 180.783, 185.133,  54.467, 143.350}
            }, { // ddg index 10: C2' C2' trans
                {149.863, 247.562, 170.488, 277.938,  84.425, 176.413, 148.087},
                {143.940, 258.200, 298.240, 279.640, 183.680, 183.080, 145.120}
            }, { // ddg index 11: C2' C2' minus
                {147.342, 256.475, 295.508, 287.408, 194.525, 293.725, 150.458}
        } };

        // Get number of ddg clusters and number of suites in each ddg cluster
        // Dominant suites lists indices of the dominant suite associated with
        // a satellite suite. A value of reference_suite_dihedrals.size()
        // that this suite is neither dominant nor satellite. A dominant suite
        // will point to its own index.
        N_reference_ddg = reference_dihedrals.size();
        N_reference_suite.clear();
        dominant_suites.clear();
        for (size_t i = 0; i < N_reference_ddg; ++i) {

            suite_size = reference_dihedrals[i].size();
            N_reference_suite.push_back(suite_size);
            vector<size_t> dom_suites(suite_size, suite_size);
            dominant_suites.push_back(dom_suites);

        }

        // 1m, 1L, and &a are satellites of 1a
        dominant_suites[0][0] = 0;
        dominant_suites[0][1] = 0;
        dominant_suites[0][2] = 0;
        dominant_suites[0][3] = 0;

        // 1f is a satellite of 1c
        dominant_suites[1][1] = 1;
        dominant_suites[1][2] = 1;

        // 1[ is a satellite of 1b
        dominant_suites[3][0] = 0;
        dominant_suites[3][1] = 0;

        // 4a and #a are satellites of 0a
        dominant_suites[6][2] = 2;
        dominant_suites[6][1] = 2;
        dominant_suites[6][3] = 2;

        // 0i nd 6j are satellites of 6n
        dominant_suites[7][3] = 3;
        dominant_suites[7][2] = 3;
        dominant_suites[7][4] = 3;

        // Two-character suite name
        reference_names = {
            {"1a", "1m", "1L", "&a", "7a", "3a", "9a", "1g", "7d", "3d", "5d",
                "3g"},
            {"1e", "1c", "1f", "5j", "5n"},
            { },
            {"1b", "1[", "3b", "1z", "5z", "7p", "5p"},
            {"1t", "5q"},
            {"1o", "7r", "5r"},
            {"2a", "4a", "0a", "#a", "4g", "6g", "8d", "4d", "6d", "2g"},
            {"2h", "4n", "0i", "6n", "6j"},
            {"0k"},
            {"2[", "4b", "0b", "4p", "6p", "2z"},
            {"4s", "2u"},
            {"2o"}
        };

        // Widths used to scale each dihedral dimension
        dihedral_width = {28.0, 60.0, 55.0, 50.0, 70.0, 35.0, 28.0};

        // Alternative widths used to scale dominant-satellite pairs
        dominant_width = {
            {28.0, 60.0, 55.0, 50.0, 64.0, 35.0, 28.0},
            {28.0, 70.0, 55.0, 50.0, 70.0, 35.0, 28.0},
            {28.0, 60.0, 60.0, 50.0, 70.0, 35.0, 28.0},
            {28.0, 60.0, 55.0, 50.0, 65.0, 35.0, 28.0},
            {28.0, 60.0, 55.0, 50.0, 56.0, 35.0, 28.0},
            {28.0, 50.0, 50.0, 50.0, 70.0, 35.0, 28.0},
            {28.0, 36.0, 36.0, 50.0, 70.0, 35.0, 28.0},
            {28.0, 60.0, 55.0, 50.0, 60.0, 35.0, 28.0},
            {28.0, 60.0, 55.0, 50.0, 60.0, 35.0, 28.0}
        };

        satellite_width = {
            {28.0, 60.0, 55.0, 50.0, 32.0, 35.0, 28.0},
            {28.0, 18.0, 55.0, 50.0, 18.0, 35.0, 28.0},
            {28.0, 20.0, 20.0, 50.0, 70.0, 35.0, 28.0},
            {28.0, 60.0, 55.0, 50.0, 47.0, 35.0, 28.0},
            {28.0, 60.0, 55.0, 50.0, 34.0, 35.0, 28.0},
            {28.0, 40.0, 40.0, 50.0, 70.0, 35.0, 28.0},
            {28.0, 26.0, 26.0, 50.0, 70.0, 35.0, 28.0},
            {28.0, 60.0, 55.0, 50.0, 60.0, 35.0, 28.0},
            {28.0, 60.0, 55.0, 50.0, 60.0, 35.0, 28.0},
        };

        // Index into dominant-satellite pair widths
        dom_sat_pair_index = {
            {9, 0, 1, 2, 9, 9, 9, 9, 9, 9, 9, 9},
            {9, 9, 3, 9, 9},
            { },
            {9, 4, 9, 9, 9, 9, 9},
            {9, 9},
            {9, 9, 9},
            {9, 5, 9, 6, 9, 9, 9, 9, 9, 9},
            {9, 9, 7, 9, 8},
            {9},
            {9, 9, 9, 9, 9, 9},
            {9, 9},
            {9}
        };

        // Delta(i-1), delta, gamma index. Delta can be C3' endo ("3") or
        // C2' endo ("2"). Gamma can be plus ("p"), trans ("t"), or minus ("m").
        reference_ddgs = {"33p", "33t", "33m", "32p", "32t", "32m", "23p",
            "23t", "23m", "22p", "22t", "22m"};

        // Boundaries for allowed regions of delta(i-1), delta, and gamma
        delta_min = { 60.0, 125.0};
        delta_max = {105.0, 165.0};
        gamma_min = { 20.0, 140.0, 260.0};
        gamma_max = { 95.0, 215.0, 335.0};

        // Boundaries used to filter suites based on epsilon, zeta, alpha, beta
        filter_min = {155.0,  25.0,  25.0,  50.0};
        filter_max = {310.0, 335.0, 335.0, 290.0};

    } // defineSuitesFromSuitename()

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

    vector<vector<double>> RnaSuite::getSuiteDihedrals() const {
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

