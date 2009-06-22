/*
  subsetter.cpp

  A general purpose tool for subsetting a trajectory.  This tool can
  be used to extract specific atoms from a trajectory or specific
  frames.  It can also be used to add periodic box information (or
  correct it) to a trajectory.  It can also be used to concatenate
  trajectories together (optionally extracting a subset of the
  concatenated trajectory).  Finally, you can center the output so
  that the centroid of the selection is the origin.  Note that the
  selection used for centering comes from the specified subset...

  The output is always in DCD format.

  Usage:
    subsetter [options] output-prefix input-model input-trajectory [input-trajectory ...]

  Examples:

  * subsetter -S10 out.dcd model.pdb traj1.dcd traj2.dcd traj3.dcd
    This concatenates the 3 trajectories together and outputs every
    10th frame to out.dcd

  * subsetter -c 'name == "CA"' out.dcd model.pdb traj1.dcd traj2.dcd traj3.dcd
    This concatenates the 3 trajectories together centering the output
    using the centroid of all c-alphas.

  * subsetter -c 'segid == "HEME"' -s '!hydrogen' out.dcd model.pdb traj.dcd
    This pulls all non-hydrogen atoms out of the trajectory and writes
    them to out.dcd, centering so that the HEME segment is at the
    origin.

  * subsetter -r 0:49 -r 150:10:300 out.dcd model.pdb traj1.dcd traj2.dcd
    This concatenates the two trajectories together, then writes out
    the first 50 frames, then frames 150 through 300 stepping by 10
    frames.  The frame indices written are of the composite
    trajectory.

  * subsetter --sort out.dcd model.pdb frames_*.dcd
    This will concatenate all frames together, sorting them
    numerically so that frames_0.dcd is first, followed by
    frames_1.dcd, frames_2.dcd, etc.

  * subsetter --sort --scanf 'run_13_%u.dcd' out.dcd model.pdb *.dcd
    This will concatenate all frames together, sorting them
    numerically as above, but will extract the second number from the
    filename as the trajectory file index.  Alternatively, the
    following option could be used in lieu of the --scanf option:
      --regex 'run_\d+_(\d+).dcd'


  Notes:

    The sorting option addresses a problem where you want to combine a
    set of trajectories that have have a linearly increasing id
    associated with them, i.e. "traj.0.dcd", "traj.1.dcd", etc.  If
    you give "traj.*.dcd" on the command-line, you will [most likely]
    get the files sorted in lexical order, not numerical order:
      traj.0.dcd
      traj.1.dcd
      traj.10.dcd
      traj.11.dcd
      ...
      traj.2.dcd
      traj.20.dcd
      ...

    Giving subsetter the "--sort" option causes subsetter to extract a
    number from the trajectory filename and sort based on that
    number.  There are two ways you can tell subsetter how to extract
    that number.  The first is to use a scanf-style format string, the
    second is to use a regular expression.  The default is to use a
    regular expression that extracts the longest sequence of digits
    from the filename...  In all cases, there is only one number that
    can be extracted and sorted on (i.e. you cannot do a two-column
    sort).

    scanf-style format

      For more detailed information, see the man-page for scanf.  In
    brief, you will want to insert a "%u" wherever the number appears
    in the filename.  In the case that you have two varying numbers,
    but you want to extract the second (or later one), use "%*u" to
    match a number without extracting it, i.e. "run_%*u_chunk_%u.dcd"

    regular expression format

      The regular expression (regex) format supported by subsetter is
    the BOOST regular expression library standard with PERL
    extensions.  The extractor looks for the first matched
    subexpression where the entire match can be converted to a
    number.  This means you can have multiple subexpressions, so long
    as the first one that is entirely a number is the one you want to
    extract one.  The default regex is "(\d+)" which means it will
    match the longest string of digits in the filename.  As in the
    example above, to match the second set of digits, use a regular
    expression like "run_\d+_(\d+).dcd".
*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
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
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <sstream>

#include <stdlib.h>

using namespace std;
using namespace loos;
namespace po = boost::program_options;

// Globals...yuck...

vector<uint> indices;            // Global indices of frames to extract
vector<uint> file_binding;       // Links global frame # to the corresponding filename
vector<uint> local_indices;      // Maps global frame numbers into local frame numbers
uint stride = 0;                 // Step through frames by this amount
                                 // (unless ranges were specified)
string model_name, out_name;
vector<string> traj_names;
string selection;
bool verbose = false;
uint verbose_updates;            // Frequency of writing updates with
                                 // verbose logging..

bool box_override = false;
GCoord box;

string center_selection;
bool center_flag = false;



// Code required for parsing trajectory filenames...

struct ScanfFmt {
  ScanfFmt(const string& s) : fmt(s) { }

  uint operator()(const string& s) const {
    uint d;
    if (sscanf(s.c_str(), fmt.c_str(), &d) != 1) {
      cerr << boost::format("Bad conversion of '%s' using format '%s'\n") % s % fmt;
      exit(-20);
    }

    return(d);
  }

  string fmt;
};



struct RegexFmt {
  RegexFmt(const string& s) : fmt(s) {
    regexp = boost::regex(s, boost::regex::perl);
  }
  
  uint operator()(const string& s) const {
    boost::smatch what;
    if (boost::regex_search(s, what, regexp)) {
      for (uint i=0; i<what.size(); ++i) {
        string submatch(what.str(i));
        char *p;
        if (submatch.empty())    // Should never happen?
          continue;

        uint val = strtoul(submatch.c_str(), &p, 10);
        if (*p == '\0')
          return(val);
      }
    }
    
    cerr << boost::format("Bad conversion of '%s' using regexp '%s'\n") % s % fmt;
    exit(-20);
  }

  string fmt;
  boost::regex regexp;
};



// Binding of trajectory name to the file # for sorting...
struct SortDatum {
  SortDatum(const string& s_, const uint n_) : s(s_), n(n_) { }
  SortDatum() : s(""), n(0) { }

  string s;
  uint n;
};


// Define comparison function for sort
bool operator<(const SortDatum& a, const SortDatum& b) {
  return(a.n < b.n);
}



// Given a vector of trajectory filenames, along with a functor for
// extracting the frame index from the filename, sorts it in numeric
// ascending order...

template<class FmtOp>
vector<string> sortNamesByFormat(vector<string>& names, const FmtOp& op) {
  uint n = names.size();
  vector<SortDatum> bound(n);
  for (uint i=0; i<n; ++i) {
    bound[i] = SortDatum(names[i], op(names[i]));
  }
  
  sort(bound.begin(), bound.end());

  vector<string> sorted(n);
  for (uint i=0; i<n; ++i)
    sorted[i] = bound[i].s;

  return(sorted);
}





void parseOptions(int argc, char *argv[]) {
  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("verbose,v", "Verbose output")
      ("updates,u", po::value<uint>(&verbose_updates)->default_value(100), "Frequency of verbose updates")
      ("selection,s", po::value<string>(&selection)->default_value("all"), "Subset selection")
      ("stride,S", po::value<uint>(), "Step through this number of frames in each trajectory")
      ("range,r", po::value< vector<string> >(), "Frames of the DCD to use (list of Octave-style ranges)")
      ("box,b", po::value<string>(), "Override any periodic box present with this one (a,b,c)")
      ("center,c", po::value<string>(), "Recenter the trajectory using this selection (of the subset)")
      ("sort", "Sort (numerically) the input DCD files.")
      ("scanf", po::value<string>(), "Sort using a scanf-style format string")
      ("regex", po::value<string>()->default_value("(\\d+)"), "Sort using a regular expression");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value< vector<string> >(&traj_names), "Trajectory filenames")
      ("out", po::value<string>(&out_name), "Output prefix");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("out", 1);
    p.add("model", 1);
    p.add("traj", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj") && vm.count("out"))) {
      cerr << "Usage- " << argv[0] << " [options] output-prefix model-name trajectory-name [trajectory-name ...]\n";
      cerr << generic;
      exit(-1);
    }

    if (vm.count("range")) {
      vector<string> ranges = vm["range"].as< vector<string> >();
      indices = parseRangeList<uint>(ranges);
    }

    if (vm.count("verbose"))
      verbose = true;

    if (vm.count("stride"))
      stride = vm["stride"].as<uint>();

    if (vm.count("box")) {
      string s = vm["box"].as<string>();
      istringstream is(s);
      if (!(is >> box)) {
        cerr << "ERROR - unable to convert " << s << ".  It must be in (a,b,c) format.\n";
        exit(-1);
      }
      box_override = 1;
    }

    if (vm.count("center")) {
      center_selection = vm["center"].as<string>();
      center_flag = true;
    }

    // Sort trajectory filenames if requested...
    if (vm.count("sort")) {
      if (vm.count("scanf")) {
        string fmt = vm["scanf"].as<string>();
        traj_names = sortNamesByFormat(traj_names, ScanfFmt(fmt));

      } else {

        string fmt = vm["regex"].as<string>();
        traj_names = sortNamesByFormat(traj_names, RegexFmt(fmt));

      }

    }

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}


uint getNumberOfFrames(const string& fname, AtomicGroup& model) {
  pTraj traj = createTrajectory(fname, model);
  return(traj->nframes());
}



// Build the mapping of global frame indices into individual files,
// and into the frame number within each file...  Also returns the
// total number of frames in the composite trajectory.

uint bindFilesToIndices(AtomicGroup& model) {
  uint total_frames = 0;

  for (uint j=0; j<traj_names.size(); ++j)  {
    uint n = getNumberOfFrames(traj_names[j], model);
    if (verbose)
      cout << boost::format("Trajectory \"%s\" has %d frames\n") % traj_names[j] % n;
    total_frames += n;
    for (uint i=0; i<n; ++i) {
      file_binding.push_back(j);
      local_indices.push_back(i);
    }
  }

  return(total_frames);
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  parseOptions(argc, argv);

  AtomicGroup model = createSystem(model_name);
  AtomicGroup subset = selectAtoms(model, selection);
  subset.clearBonds();

  AtomicGroup centered;
  if (!center_selection.empty())
    centered = selectAtoms(subset, center_selection);

  uint total_frames = bindFilesToIndices(model);

  // If no frame ranges were specified, fill in all possible frames,
  // using the stride (if given)...
  if (indices.empty()) {
    stride = (!stride) ? 1 : stride;
    for (uint i=0; i<total_frames; i += stride)
      indices.push_back(i);
  }

  DCDWriter dcdout(out_name + ".dcd");
  dcdout.setTitle(hdr);

  bool first = true;  // Flag to pick off the first frame for a
                      // reference structure
  pTraj traj;
  uint current = 0;   // Track the index of the current trajectory
                      // we're reading from...

  // Setup for progress output...
  PercentProgressWithTime watcher;
  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(PercentTrigger(0.25), EstimatingCounter(indices.size()));
  slayer.attach(&watcher);
  if (verbose)
    slayer.start();

  // Iterate over all requested global-frames...
  vector<uint>::iterator vi;
  for (vi = indices.begin(); vi != indices.end(); ++vi) {

    // Have we switched to a new file??
    // Note: Because of the way file bindings are setup, it is
    // possible to get an out-of-range index prior to trying to read
    // the trajectory frame...  This is solved by using the at()
    // method here.
    if (file_binding.at(*vi) != current || first) {
      current = file_binding[*vi];
      traj = createTrajectory(traj_names[current], model);
    }

    // Read the apropriate local frame...
    traj->readFrame(local_indices[*vi]);
    traj->updateGroupCoords(model);

    // Handle centering...
    if (center_flag) {
      GCoord c = centered.centroid();
      model.translate(-c);
    }

    // Handle Periodic boundary conditions...
    if (box_override) {
      if (first && subset.isPeriodic())
        cerr << "WARNING - overriding existing periodic box.\n";
      subset.periodicBox(box);
    }

    dcdout.writeFrame(subset);

    // Pick off the first frame for the reference structure...
    if (first) {
      PDB pdb = PDB::fromAtomicGroup(subset);
      if (selection != "all")     // Strip connectivity if subsetting,
        pdb.clearBonds();         // otherwise the PDB writer will fail...
      pdb.remarks().add(hdr);
      string out_pdb_name = out_name + ".pdb";
      ofstream ofs(out_pdb_name.c_str());
      ofs << pdb;
      ofs.close();
      first = false;
    }

    if (verbose)
      slayer.update();
  }

  if (verbose)
    slayer.finish();
}
