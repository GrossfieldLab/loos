/*
  long-bond-finder.cpp

  (c) 2021 Louis Smith, Bowman Lab
           Department of Biochemistry & Biophysics
           Washington University in St. Louis, School of Medicine


  A tool that identifies trajectories, and optionally frames/atoms,
  that are overlong--a way to find imaging issues and distorted structures.
*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011 Tod D. Romo
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

#include <iostream>
#include <loos.hpp>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

// clang-format off
const string msg = 
"SYNOPSIS \n"
" \n"
"This tool is designed to find trajectories that have distorted structures or \n"
"bad rewrapping by checking for overlong bonds in some frame within the \n"
"trajectory. There are two modes of operation. In the first, the tool scans all \n"
"frames until it finds a problem, whereupon it either returns 1 (EXIT_FAILURE) \n"
"if it has, or 0 (EXIT_SUCCESS) if not. The second additionally writes a time \n"
"series to a provided file name listing each violating bond for each frame.  \n"
"Note that merge-traj with --fix-imaging may be able to resolve issues flagged\n"
"by this program. gmx-trjconv with cluster imaging may also work.\n"
" \n"
"DESCRIPTION \n"
" \n"
"The tool uses predefined connectivity to look through all frames for bonds that\n"
" are longer than some cutoff within the user-specified selection. If it is \n"
"operating in default mode, that is if no time-series is requested, it will exit\n"
" upon finding the first such flawed bond. It will either return an \n"
"'EXIT_FAILURE' if it has found at least one, or 'EXIT_SUCCESS' if not. It can \n"
"also write a time-series of these bonds to a filename provided by the user, in \n"
"which case it will check all frames and report any flawed bonds to the file \n"
"name provided. For processing large datasets the first mode is likely more \n"
"helpful and will almost certainly be faster. \n"
" \n"
"The intention is for the return value to be used in bash control flow to allow \n"
"a user to conditionally add a trajectory file to a list of files with problems \n"
"(probably PBC issues, but I don't know your life), or perhaps operate on them \n"
"directly if the tool flags them (see 'EXAMPLES' below). This is analogous to \n"
"using 'grep -q' (quiet mode) to ascertain whether a regex is somewhere inside a\n"
" file,  then operating conditionally in response to whether it is present or \n"
"absent. The tool provides a quiet mode to support this scripting style. The \n"
"time-series mode of operation can be nice as a way to spot check where the \n"
"issue is, since any frames noted there can then be looked up in a visualizer \n"
"and inspected manually to see what is really going wrong. The the return-value \n"
"of the program works the same whether a time-series is requested or not. \n"
" \n"
"If a time-series filename is provided as an argument to the option, a comment-\n"
"line containing the invocation is written to the first line.  Subsequent lines \n"
"have four values separated by a space, each representing a bond that is \n"
"overlong: the frame index, the first atomID within the pair of atoms \n"
"constituting some bond, the second atomID in that bond, and that bond's \n"
"calculated length. Only non-redundant ID pairs in an order-independent fashion \n"
"are checked or reported (since atomID 2 bonded to atomID 7 is the same as 7 \n"
"bonded to 2).  \n"
" \n"
"EXAMPLES \n"
" \n"
"The most basic mode of operation for this tool is: \n"
" \n"
"long-bond-finder model.psf traj.dcd \n"
"if [[ $? > 0 ]]; then \n"
"  echo this traj is goofed \n"
"fi \n"
" \n"
"More useful is the situation where this is used in some loop, silencing the \n"
"per-traj operations using output redirection: \n"
" \n"
"for traj in xtcs/*.xtc; do \n"
"  if ! long-bond-finder model.psf $traj 1 >> logfile.log 2>&1; then \n"
"    echo $traj is goofed. 'DO SOMETHING!!!' \n"
"    echo $traj >> list_of_goofy_trajs.txt \n"
"    mv $traj goofy_trajs/ \n"
"  fi \n"
"done \n"
" \n"
"If no log-file is desired, throw the '--quiet' flag to suppress emission of the\n"
" invocation header. \n"
" \n"
"for traj in xtcs/*.xtc; do \n"
"  if ! long-bond-finder --quiet model.psf $traj; then \n"
"    echo $traj is goofed. 'DO SOMETHING!!!' \n"
"    echo $traj >> list_of_goofy_trajs.txt \n"
"    mv $traj goofy_trajs/ \n"
"  fi \n"
"done \n"
" \n"
" \n"
"Running either of these two commands with the addition of the '--timeseries \n"
"filename.dat' flag will write the four-column timeseries of bond violations \n"
"(possibly none) to 'filename.dat'. \n"
" \n"
"To infer connectivity for some model, provide the cutoff distance for a 'bond' \n"
"as an argument to the option: \n"
" \n"
"long-bond-finder --infer-connectivity 1.9 my_minimal_model.pdb traj.dcd \n"
" \n"
"To change what the cutoff for distorted bonds is, use the '--max-bond value' \n"
"flag. For example, if a 3.2 Angstrom cutoff were desired: \n"
" \n"
"long-bond-finder --max-bond 3.2 model.psf traj.dcd \n"
" \n"
"POTENTIAL COMPLICATIONS  \n"
" \n"
"Note that as with all loos tools, trajectory file indexes are zeros based, but \n"
"some visualizers can be ones-based. In addition, GROMACS-generated trajectories\n"
" sometimes have initial coords saved to the 'first' frame, leaving the \n"
"possibility that the literal value in the first column could be off by either \n"
"one or two from what is displayed in a visualizer, depending on the \n"
"circumstances. For issue diagnosis purposes this seems OK, and the zeros-based \n"
"index is correct if one wanted to write a secondary script that used those \n"
"values as indexes into the trajectory to do something to goofed up frames. \n"
" \n"
"This tool doesn't currently use periodicity for bond length calculations, in \n"
"part because it was written to spot bad wrapping issues from harder to wrap \n"
"systems with periodicity that is not presently supported in loos (rhombic \n"
"dodecahedra, for example). Because many loos tools do use periodicity for \n"
"distance calculations, this could surprise some users. It seems hard to catch \n"
"bad wrapping with a distance calculation if that calculation respected PBCs. \n"
" \n"
"Also note that although a mechanism is provided to use models that don't have \n"
"connectivity, this option should be deployed cautiously. It uses a simple \n"
"distance cutoff to deduce where chemical bonds likely are for the model as a \n"
"whole based on the coordinates in the first frame. This will be incorrect if \n"
"the first frame has bonds that are overlong (or extreme collisions) relative to\n"
" the user's expectation, and while the tool is likely to find something \n"
"objectionable in the successive screwed up mess that its output will become, \n"
"manual inspection would be needed to be sure that the objection was not a false\n"
" positive. Using models with chemical connectivity based on a more reliable \n"
"source, such as a system-specifying file from an MD engine, is highly advised. \n"
" \n"
"Finally, note that the time-series mode could produce unwieldy output in the \n"
"case where a system specifying file is fundamentally damaged (say the thing \n"
"that is wrong with it is that its bond indices are off by one for all atoms). \n"
"The potential for absurdly voluminous output can be guarded against by checking\n"
" whether a trajectory has major issues by not writing the time-series as a \n"
"first past, then if a trajectory is flagged visualizing its contents to ensure \n"
"that the time-series won't have a ridiculous number (possibly nearly all) of \n"
"bonds to record for each frame. Similarly, using the customary loos trajectory \n"
"flags (range) alongside the time-series flag could allow one to only consider \n"
"the report on the first few frames, then inspecting these for pathology which \n"
"would make seeing the rest of the time-series irrelevant and overwhelming. \n"
;
// clang-format on

// The following conditional prevents this class from appearing in the
// Doxygen documentation for LOOS:
//
// @cond TOOLS_INTERNAL
// clang-format off
class ToolOptions : public opts::OptionsPackage {
public:

  // Change these options to reflect what your tool needs
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("infer-connectivity", po::value<float>(&bondlength)->default_value(-1), 
       "Infer connectivity using provided distance for models lacking this. "
       "ALERT: uses provided value as hard distance cutoff on first frame of traj to infer connectivity. "
       "Only does this for values greater than zero.")
      ("max-bond,M", po::value<float>(&max_bond)->default_value(2.5),
       "Maximum permissible distance for plausible bond.")
      ("quiet,q", po::bool_switch(&quiet)->default_value(false),
       "Silence standard output.")
      ("timeseries,t", po::value<string>(&timeseries)->default_value(""), 
       "Write bond-pairs in violation of cutoff per-frame to file name provided.");
  }

  // The print() function returns a string that describes what all the
  // options are set to (for logging purposes)
  string print() const {
    ostringstream oss;
    oss << boost::format("bondlength=%f,max_bond=%f,quiet=$%b,timeseries=%s") % bondlength % max_bond % quiet % timeseries;
    return(oss.str());
  }

  float bondlength, max_bond;
  bool quiet;
  string timeseries;

};
// clang-format on
// @endcond
// ----------------------------------------------------------------

int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);
  opts::BasicOptions *bopts = new opts::BasicOptions(msg);
  opts::BasicSelection *sopts = new opts::BasicSelection("all");
  opts::TrajectoryWithFrameIndices *tropts =
      new opts::TrajectoryWithFrameIndices;
  ToolOptions *topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts);

  if (!options.parse(argc, argv))
    exit(-1);

  // set up system for looping. Load coords from frame 0 into scope.
  AtomicGroup model = tropts->model;
  AtomicGroup scope = selectAtoms(model, sopts->selection);
  bool all_bonds_in_scope = scope.allHaveProperty(Atom::bits::bondsbit);
  if (all_bonds_in_scope) {
  } else if (topts->bondlength > 0)
    if (scope.hasCoords())
      scope.findBonds(topts->bondlength);
    else {
      throw(LOOSError(
        "Model does not have coordinates with which to infer connectivity.\n"
      ));
    }
  else
    throw(LOOSError(
        "Model selection does not appear to have chemical connectivity, and "
        "infer-connectivity has not been set to a positive value.\n"));
  pTraj traj = tropts->trajectory;
  traj->updateGroupCoords(model);
  // should be a vector of two-atom AGs, each a pair of atoms in a bond
  vector<AtomicGroup> bond_list = scope.getBondsAGs();
  // set threshold for length of an unacceptable bond.
  const double max_bond2 = topts->max_bond * topts->max_bond;

  // Operating in scanning mode;
  // don't report anything except the presence of an unacceptable bond
  if (topts->timeseries.empty()) {
    // if thrown, don't even write invocation to stdout
    if (!topts->quiet)
      cout << "# " << header << "\n";
    
    float dist2 = 0;
    for (auto frame_index : tropts->frameList()) {
      traj->readFrame(frame_index);
      traj->updateGroupCoords(scope);
      for (auto b : bond_list) {
        dist2 = b[0]->coords().distance2(b[1]->coords());
        if (dist2 > max_bond2) {
          if (!topts->quiet)
            cout << "Issue in frame " << frame_index << "; bond between atomIDs " << b[0]->id() << " and " << b[1]->id() << 
            " is " << sqrtf(dist2) << " Angstroms. Exiting..." << endl;
          return EXIT_FAILURE;
        }
      }
    }
  } else { 
    // Operating in timeseries mode;
    // write a timeseries to file name provided by user
    ofstream tsf(topts->timeseries);
    bool found_viol = false;
    tsf << "# " << header << "\n"
        << "# frame atomID1 atomID2 bondlength\n";
    for (auto frame_index : tropts->frameList()) {
      traj->readFrame(frame_index);
      traj->updateGroupCoords(scope);
      float dist2 = 0;
      for (auto b : bond_list) {
        dist2 = b[0]->coords().distance2(b[1]->coords());
        if (dist2 > max_bond2) {
          found_viol = true;
          tsf << frame_index << " " << b[0]->id() << " " << b[1]->id() 
              << " " << sqrtf(dist2) << "\n";
        }
      }
    }
    if (found_viol)
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
