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
"\n"
"    SYNOPSIS\n"
"\n"
"    Assigns backbone suites to RNAs based on backbone dihedrals\n"
"\n"
"    DESCRIPTION\n"
"\n"
"    This tool\n"
"\n"
"    EXAMPLES\n"
"\n"
"    rna_suites\n"
        ;

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
    opts::TrajectoryWithFrameIndices *tropts =
        new opts::TrajectoryWithFrameIndices;
    ToolOptions *topts = new ToolOptions;

    opts::AggregateOptions options;
    options.add(bopts).add(sopts).add(tropts).add(topts);
    if (!options.parse(argc, argv))
        exit(-1);

    // Assign tool options to variables
    const double suiteness_cutoff = topts->suiteness_cutoff;

    // Print command-line input
    cout << "# " << header << "\n";

    // Do some error-checking on tool options

    // Build LOOS system and generate atom selection
    AtomicGroup model = tropts->model;
    pTraj traj = tropts->trajectory;
    vector<uint> indices = tropts->frameList();
//    AtomicGroup rna_atoms = selectAtoms(model, topts->selection);

    // Number of frames in trajectory
    const uint N_frame = indices.size();

    // Create RNASuite object from RNA atoms
    RnaSuite rna_suite = RnaSuite(model, suiteness_cutoff);

    // Print dihedrals
    rna_suite.printBackboneAtoms();

}
