/*
  rna_suites.cpp

  Assigns backbone suites to RNAs based on backbone dihedrals

  Chapin E. Cavender 2020-03
*/

/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008-2020 Tod D. Romo & Alan Grossfield
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

#include <loos.hpp>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

string fullHelpMessage(void) {

    string full_help_message =
"SYNOPSIS\n"
"    Assign backbone suites to RNAs based on backbone dihedrals.\n"
"\n"
"DESCRIPTION\n"
"    The goal of this tool is to assign continuous RNA dinucleotides to a\n"
"    cluster called a \"suite\" based on the conformation of backbone dihedrals.\n"
"    The idea comes from Richardson et al. (2008) RNA 14, 465-481. The\n"
"    dinucleotide for a residue runs from delta (C5'-C4'-C3'-O3') of the previous\n"
"    residue to delta of the current residue, encompassing seven continuous\n"
"    dihedrals. A suite is a pre-defined cluster in this 7D space, named by a\n"
"    two-character string. Examples are \"1a\" or \"5z\".\n"
"\n"
"    The first step is to search the given selection for RNA backbone atoms, i.e.\n"
"    atoms named \"P\", \"O5'\", \"C5'\", \"C4'\", \"C3'\", or \"O3'\". These atoms are\n"
"    split by residue. Valid dinucleotides are sets of delta-to-delta backbone\n"
"    atoms with sequential resids. Once the set of valid dinucleotides is\n"
"    determined, the tool will loop over the trajectory and assign each\n"
"    dinucleotide to a suite for each frame.\n"
"\n"
"    Suite assignment occurs in two stages. The clusters are well-separated in\n"
"    the 3D subspace of delta(i-1), delta, and gamma. So the first stage is to\n"
"    assign each delta to one of two ranges of values consistent with either a\n"
"    C3'-endo (3) or C2'-endo (2) sugar pucker and to assign gamma to one of\n"
"    three ranges of values: gauche plus (p), gauche minus (m), or trans (t). The\n"
"    result is a three-character string called a ddg index. Examples are \"33p\"\n"
"    or \"23t\". Then, the dinucleotide is assigned to one of a possible set of\n"
"    suites associated with its ddg index based on a scaled hyperellipsoid\n"
"    distance in the dual 4D subspace of epsilon, zeta, alpha, and beta.\n"
"\n"
"    Some suites have overlapping hyperellipsoids of different sizes. The wider\n"
"    suite is called a dominant suite, and the narrower suite is called a\n"
"    satellite suite. These cases are handled by rescaling the hyperellipsoid\n"
"    distance along the dimensions in which the overlap occurs.\n"
"\n"
"    If a dinucleotide doesn't fit into one of the allowed ranges for a dihedral,\n"
"    it is assigned as an outlier and given a suite name \"!s\", where \"s\" is the\n"
"    first character of the name of the deviant dihedral, e.g. \"!a\" for a bad\n"
"    alpha. If the dinucleotide is not close to any of the reference suites, it \n"
"    is also assigned as an outlier and given a suite name \"!!\".\n"
"\n"
"    After assignment, each dinucleotide is given a goodness-of-fit score called \n"
"    the suiteness based on the scaled 7D hyperellipsoid distance to its assigned\n"
"    suite. A suiteness of one indicates that the dinucleotide is at the cluster\n"
"    center. Lower suiteness indicates that the dinucleotide is farther from the\n"
"    cluster center. An outlier has a suiteness of zero, and assigned\n"
"    dinucleotides have a minimum suiteness score (set by the -c option) to\n"
"    differentiate them from outliers.\n"
"\n"
"    It is necessary to specify a path to a file containing definitions for the\n"
"    reference suites on the command-line. The format is explained in the next\n"
"    section. An example of the format that implements the suites as defined in\n"
"    the software suitename (Richardson et al. (2008) RNA 14, 465-481) is\n"
"    included as share/suitename_definitions.dat in the top-level directory of\n"
"    the LOOS source tree. If installing within a conda environment, this file\n"
"    can also be found in $CONDA_PREFIX/share/loos/suitename_definitions.dat;\n"
"    otherwise, it can be found in $LOOS/share/suitename_definitions.dat. The\n"
"    suitename_defintions.dat file should be sufficient for typical users, but\n"
"    you must specify the path to it as the first positional argument.\n"
"\n"
"SUITE DEFINITON FILE FORMAT\n"
"    Each line in the file is parsed as a record containing fields with a width\n"
"    of eight characters. Blank lines and lines beginning with \"#\" are ignored.\n"
"    The first field specifies the type of record and must be one of \"suite\",\n"
"    \"width\", \"domsat\", \"delta\", \"epsilon\", \"zeta\", \"alpha\", \"beta\", or \"gamma\".\n"
"    These records and their associated fields are described below.\n"
"\n"
"    suite name ddg delta(i-1) epsilon zeta alpha beta gamma delta(i)\n"
"        Define a reference suite with suite name given in field 2, ddg index\n"
"        given in field 3, and dihedrals of the cluster center given in fields 4\n"
"        through 10.\n"
"\n"
"    width delta(i-1) epsilon zeta alpha beta gamma delta\n"
"        Define default widths for scaled hyperellipsoid distances.\n"
"\n"
"    domsat sat_name dom_name dihedral_index sat_width dom_width\n"
"        Define dominant-satellite pair with name of satellite suite in field 2,\n"
"        name of dominant suite in field 3, index of dihedral dimension with\n"
"        altered width in field 4, width of that dimension for satellite suite\n"
"        in field 5, and width of that dimension for dominant suite in field 6.\n"
"        Additional dimensions and widths can be specified in fields 7 through 9,\n"
"        fields 10 through 12, etc.\n"
"\n"
"    dihedral min max\n"
"        Define allowed ranges for a dihedral. \"dihedral\" can be one of \"delta\",\n"
"        \"epsilon\", \"zeta\", \"alpha\", \"beta\", or \"gamma\". The minimum value\n"
"        is given in field 2 and maximum value in field 3.\n"
"\n"
"EXAMPLES\n"
"    rna_suites $CONDA_PREFIX/share/loos/suitename_defintions.dat foo.pdb foo.dcd\n"
"        Assign backbone suites using the install prefix from a conda install.\n"
"\n"
"    rna_suites -s 'resid <= 10' $CONDA_PREFIX/share/loos/suitename_defintions.dat \\\n"
"        foo.pdb foo.dcd\n"
"        Assign backbone suites only for the first 10 residues.\n"
"\n"
"    rna_suites -c 0.001 $CONDA_PREFIX/share/loos/suitename_defintions.dat \\\n"
"        foo.pdb foo.dcd\n"
"        Assign backbone suites using a minimum suiteness of 0.001 for\n"
"        non-outliers.\n";

    return full_help_message;

}

class ToolOptions : public opts::OptionsPackage {

public:

    ToolOptions() {}

    void addGeneric(po::options_description& o) {

        o.add_options()
            ("suiteness_cutoff,c",
                po::value<double>(&suiteness_cutoff)->default_value(0.01),
                "Cutoff for the suiteness score of non-outliers")
            ;

    }

    string print() const {

        ostringstream oss;
        oss << boost::format(
            "suiteness_cutoff=%f"
            ) % suiteness_cutoff;
        return (oss.str());

    }

    double suiteness_cutoff;

}; // ToolOptions

// Tool functions

int main(int argc, char *argv[]) {

    // Get command-line input
    string header = invocationHeader(argc, argv);

    // Set up tool options
    opts::BasicOptions *bopts = new opts::BasicOptions(fullHelpMessage());
    opts::BasicSelection *sopts = new opts::BasicSelection("!hydrogen");
    opts::RequiredArguments *ropts = new opts::RequiredArguments;
    ropts->addArgument("suite_def", "suite_definition_file");
    opts::TrajectoryWithFrameIndices *tropts =
        new opts::TrajectoryWithFrameIndices;
    ToolOptions *topts = new ToolOptions;

    opts::AggregateOptions options;
    options.add(bopts).add(sopts).add(ropts).add(tropts).add(topts);
    if (!options.parse(argc, argv))
        exit(-1);

    // Assign tool options to variables
    const double suiteness_cutoff = topts->suiteness_cutoff;

    // Print command-line input
    cout << "# " << header << endl;

    // Build LOOS system and generate atom selection
    AtomicGroup model = tropts->model;
    pTraj traj = tropts->trajectory;
    vector<uint> indices = tropts->frameList();
    AtomicGroup rna_atoms = selectAtoms(model, sopts->selection);

    // Create RNASuite object from RNA atoms
    string suite_definition = ropts->value("suite_def");
    RnaSuite rna_suite = RnaSuite(rna_atoms, suite_definition, suiteness_cutoff);
    vector<int> suite_resids = rna_suite.getSuiteResids();
    vector<string> suite_resnames = rna_suite.getSuiteResnames();
    //rna_suite.printReferenceSuites();

    // Print dihedrals
    //rna_suite.printBackboneAtoms();

    // Print column headers
    cout << "# Frame Resid Resname Suite DDG_index Suiteness" << endl;

    // Loop over trajectory
    vector<string> suite_names;
    vector<string> suite_ddgs;
    vector<double> suiteness;
    uint t = 0;
    for (vector<uint>::iterator i = indices.begin(); i != indices.end(); ++i) {

        traj->readFrame(*i);
        traj->updateGroupCoords(model);

        rna_suite.calculateBackboneDihedrals();
        rna_suite.assignSuitenameSuites();
        suite_names = rna_suite.getSuiteNames();
        suite_ddgs = rna_suite.getSuiteDDGs();
        suiteness = rna_suite.getSuitenessScores();

        for (uint j = 0; j < suite_resids.size(); ++j)
            cout << boost::format("%5d %5d %3s %2s %2s %8.6f") % t
                % suite_resids[j] % suite_resnames[j] % suite_names[j]
                % suite_ddgs[j] % suiteness[j] << endl;

        ++t;

    }

}

