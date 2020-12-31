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
#include <boost/regex.hpp>
#include <boost/lambda/lambda.hpp>
#include <sstream>

#include <cstdlib>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


// @cond TOOLS_INTERNAL


typedef vector<AtomicGroup>   vGroup;


// Globals...yuck...

vector<uint> indices;            // Global indices of frames to extract
vector<uint> file_binding;       // Links global frame # to the corresponding filename
vector<uint> local_indices;      // Maps global frame numbers into local frame numbers
uint stride = 0;                 // Step through frames by this amount
                                 // (unless ranges were specified)
uint skip = 0;
string model_name, out_name;
vector<string> traj_names;
string selection;
int verbose = 0;
uint verbose_updates;            // Frequency of writing updates with
                                 // verbose logging..

bool box_override = false;
GCoord box;


enum ReimageMode { NONE, NORMAL, AGGRESSIVE, ZEALOUS, EXTREME } reimage_mode;

ulong extreme_iters = 0;
double extreme_delta = 0.0;
const uint extreme_max_iters = 250;
const double extreme_threshold = 1e-1;

string center_selection;
bool center_flag = false;
string post_center_selection;



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





string fullHelpMessage(void) {

  string s =
    "\n"
    "SYNOPSIS\n"
    "\tConversion of trajectories to DCD format and extraction of subsets\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis is a general-purpose tool (similar to catdcd from NAMD) that can be used\n"
    "to convert trajectories to the DCD format, extract ranges of frames from a trajectory,\n"
    "extract only a subset of atoms, assign a periodic box, reimage (for periodic boundaries),\n"
    "and center the system, among others.\n"
    "\n"
    "\tReimaging can be handled several different ways.  The simplest is to turn on reimaging\n"
    "with --reimage=normal.  This reimages by molecule.  In some cases, this is insufficient to\n"
    "reimage the system so that all molecules are 'together'.  The second method is invoked with\n"
    "--reimage=aggressive.  This employs a more aggressive reimaging that attempts to keep all parts of\n"
    "a molecule together (the method used is similar to the --fix-imaging option).  A similar\n"
    "reimaging strategy is to use --reimage=zealous, a two-pass strategy where first normal\n"
    "reimaging is applied, followed by aggressive.  This can be helpful with some split GROMACS systems.\n"
    "An even more aggressive method, is to use --reimage=extreme.  Here, an iterative reimaging\n"
    "procedure is used.  This may slow down subsetter.  In aggressive, zealous, and extreme, a centering\n"
    "selection is used.  For aggressive and zealous, you should center on whatever you want the system to\n"
    "be centered on (e.g. a protein or a membrane).  Extreme can work with 'all' as a selection.  If\n"
    "that fails, try selecting either a central protein or a membrane.  Since these reimaging methods can\n"
    "affect the centering, a post-reimaging centering is available using the --postcenter option.\n"
    "Finally, these imaging methods require connectivity and, in the case of extreme, masses are\n"
    "helpful.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tsubsetter -S10 out model.pdb traj1.dcd traj2.dcd traj3.dcd\n"
    "This concatenates the 3 trajectories together and outputs every\n"
    "10th frame to out.dcd\n"
    "\n"
    "\tsubsetter -c 'name == \"CA\"' out model.pdb traj1.dcd traj2.dcd traj3.dcd\n"
    "This concatenates the 3 trajectories together centering the output\n"
    "using the centroid of all c-alphas.\n"
    "\n"
    "\tsubsetter -c 'segid == \"HEME\"' -s '!hydrogen' out model.pdb traj.dcd\n"
    "This pulls all non-hydrogen atoms out of the trajectory and writes\n"
    "them to out.dcd, centering so that the HEME segment is at the\n"
    "origin.\n"
    "\n"
    "\tsubsetter -r 0:49,150:10:300 out model.pdb traj1.dcd traj2.dcd\n"
    "This concatenates the two trajectories together, then writes out\n"
    "the first 50 frames, then frames 150 through 300 stepping by 10\n"
    "frames.  The frame indices written are of the composite\n"
    "trajectory.\n"
    "\n"
    "\tsubsetter --sort out model.pdb frames_*.dcd\n"
    "This will concatenate all frames together, sorting them\n"
    "numerically so that frames_0.dcd is first, followed by\n"
    "frames_1.dcd, frames_2.dcd, etc.\n"
    "\n"
    "\tsubsetter --sort --scanf 'run_13_%u.dcd' out model.pdb *.dcd\n"
    "This will concatenate all frames together, sorting them\n"
    "numerically as above, but will extract the second number from the\n"
    "filename as the trajectory file index.  Alternatively, the\n"
    "following option could be used in lieu of the --scanf option:\n"
    " --regex 'run_\\d+_(\\d+).dcd'\n"
    "\n"
    "\tsubsetter -t xtc out model.pdb *.dcd\n"
    "Writes out an XTC formatted trajectory to out.xtc and model to\n"
    "out.pdb.  Concatenates all DCD trajectories in the current\n"
    "directory."
    "\n"
    "\tsubsetter --reimage=extreme --center='all' --postcenter='segid == \"POPC\" out.dcd model.psf *.dcd\n"
    "Writes out a DCD reimaging the system using the extreme method and centering\n"
    "(after reimaging) on the POPC membrane\n"
    "\n"
    "NOTES\n"
    "\n"
    "\t* sorting *\n"
    "\tThe sorting option addresses a problem where you want to combine a\n"
    "set of trajectories that have have a linearly increasing id\n"
    "associated with them, i.e. \"traj.0.dcd\", \"traj.1.dcd\", etc.  If\n"
    "you give \"traj.*.dcd\" on the command-line, you will [most likely]\n"
    "get the files sorted in lexical order, not numerical order:\n"
    "  traj.0.dcd\n"
    "  traj.1.dcd\n"
    "  traj.10.dcd\n"
    "  traj.11.dcd\n"
    "  ...\n"
    "  traj.2.dcd\n"
    "  traj.20.dcd\n"
    "  ...\n"
    "\n"
    "\tGiving subsetter the \"--sort\" option causes subsetter to extract a\n"
    "number from the trajectory filename and sort based on that\n"
    "number.  There are two ways you can tell subsetter how to extract\n"
    "that number.  The first is to use a scanf-style format string, the\n"
    "second is to use a regular expression.  The default is to use a\n"
    "regular expression that extracts the longest sequence of digits\n"
    "from the filename...  In all cases, there is only one number that\n"
    "can be extracted and sorted on (i.e. you cannot do a two-column\n"
    "sort).\n"
    "\n"
    "\t* scanf-style format *\n"
    "For more detailed information, see the man-page for scanf.  In\n"
    "brief, you will want to insert a \"%u\" wherever the number appears\n"
    "in the filename.  In the case that you have two varying numbers,\n"
    "but you want to extract the second (or later one), use \"%*u\" to\n"
    "match a number without extracting it, i.e. \"run_%*u_chunk_%u.dcd\"\n"
    "\n"
    "\t* regular expression format\n"
    "The regular expression (regex) format supported by subsetter is\n"
    "the BOOST regular expression library standard with PERL\n"
    "extensions.  The extractor looks for the first matched\n"
    "subexpression where the entire match can be converted to a\n"
    "number.  This means you can have multiple subexpressions, so long\n"
    "as the first one that is entirely a number is the one you want to\n"
    "extract one.  The default regex is \"(\\d+)\" which means it will\n"
    "match the longest string of digits in the filename.  As in the\n"
    "example above, to match the second set of digits, use a regular\n"
    "expression like \"run_\\d+_(\\d+).dcd\".\n"
    "\n"
    "SEE ALSO\n"
    "\tmerge-traj, reimage-by-molecule, recenter-trj\n"
    "\n";

  return(s);
}


// Note: We do not use the TrajectoryWithFrameIndices class here
// because this tool supports a more complex arrangement of
// trajectories with ranges and skips...
class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("updates", po::value<uint>(&verbose_updates)->default_value(100), "Frequency of verbose updates")
      ("stride,i", po::value<uint>(&stride)->default_value(1), "Step through this number of frames in each trajectory")
      ("skip,k", po::value<uint>(&skip)->default_value(0), "Skip these frames at start of each trajectory")
      ("range,r", po::value<string>(&range_spec)->default_value(""), "Frames of the DCD to use (list of Octave-style ranges)")
      ("box,B", po::value<string>(&box_spec), "Override any periodic box present with this one (a,b,c)")
      ("reimage", po::value<string>(&reimage)->default_value("none"), "Reimage mode (none, normal, aggressive, zealous, extreme)")
      ("center,C", po::value<string>(&center_selection)->default_value(""), "Recenter the trajectory using this selection (of the subset)")
      ("postcenter,P", po::value<string>(&post_center_selection)->default_value(""), "Recenter using this selection after reimaging")
      ("sort", po::value<bool>(&sort_flag)->default_value(false), "Sort (numerically) the input DCD files.")
      ("scanf", po::value<string>(&scanf_spec)->default_value(""), "Sort using a scanf-style format string")
      ("regex", po::value<string>(&regex_spec)->default_value("(\\d+)\\D*$"), "Sort using a regular expression");
  }

  void addHidden(po::options_description& o) {
    o.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value< vector<string> >(&traj_names), "Trajectory filenames")
      ("out", po::value<string>(&out_name), "Output prefix");
  }

  void addPositional(po::positional_options_description& o) {
    o.add("out", 1);
    o.add("model", 1);
    o.add("traj", -1);
  }

  bool check(po::variables_map& vm) {
    return ( model_name.empty() || out_name.empty() || traj_names.empty());
  }

  bool postConditions(po::variables_map& vm) {
    if (!box_spec.empty()) {
      istringstream is(box_spec);
      try {
          is >> box;
      }
      catch (const std::exception& e) {
          cerr << e.what() << std::endl;
          cerr << "ERROR: unable to convert " << box_spec << ".  It must be in '(a,b,c)' format.\n";
          return(false);
      }
      box_override = 1;
    }

    if ( !(scanf_spec.empty() || regex_spec.empty() )) {
        sort_flag = true;
    }

    if (sort_flag) {
      if (!scanf_spec.empty()) {
        traj_names = sortNamesByFormat(traj_names, ScanfFmt(scanf_spec));
      } else {
        traj_names = sortNamesByFormat(traj_names, RegexFmt(regex_spec));
      }
    }

    center_flag = !center_selection.empty();


    if (boost::iequals(reimage, "none"))
      reimage_mode = NONE;
    else if (boost::iequals(reimage, "normal"))
      reimage_mode = NORMAL;
    else if (boost::iequals(reimage, "aggressive"))
      reimage_mode = AGGRESSIVE;
    else if (boost::iequals(reimage, "zealous"))
      reimage_mode = ZEALOUS;
    else if (boost::iequals(reimage, "extreme"))
      reimage_mode = EXTREME;
    else {
      cerr << "Error- '" << reimage << "' is an unknown reimaging mode.\n";
      cerr << "       Must be: none, normal, aggressive, extreme.\n";
      return(false);
    }

    if ((reimage_mode == AGGRESSIVE || reimage_mode == ZEALOUS || reimage_mode == EXTREME) && !center_flag) {
      cerr << "Error- aggressive, zealous, and extreme reimaging modes require a centering selection.\n";
      return(false);
    }

    return(true);
  }

  string help() const {
    return("output-prefix model trajectory [trajectory ...]");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("updates=%d, stride=%s, skip=%d, range='%s', box='%s', reimage='%s', center='%s', sort=%d, postcenter='%s'")
      % verbose_updates
      % stride
      % skip
      % range_spec
      % box_spec
      % reimage
      % center_selection
      % sort_flag
      % post_center_selection;
    if (sort_flag) {
      if (!scanf_spec.empty())
        oss << boost::format("scanf='%s'") % scanf_spec;
      else
        oss << boost::format("regex='%s'") % regex_spec;
    }

    oss << boost::format("out='%s', model='%s', traj='%s'")
      % out_name
      % model_name
      % vectorAsStringWithCommas(traj_names);

    return(oss.str());
  }


  bool sort_flag;
  string regex_spec;
  string scanf_spec;
  string range_spec;
  string box_spec;
  string reimage;
};




void showTrajectoryTable(MultiTrajectory& traj) {

  cout << "Input Trajectory Table:\n";

  cout << boost::format("%5s %8s %8s %8s %8s %s\n")
    % "Traj" % "Start" % "End" % "N" % "Orig" % "Name";

  cout << boost::format("%5s %8s %8s %8s %8s %s\n")
    % "----" % "-----" % "---" % "-" % "----" % "----";

  uint start_cnt = 0;
  uint j = 0;
  for (uint i=0; i<traj.size(); ++i) {
    uint n = traj.nframes(i);
    if (n == 0)
      cout << boost::format("%5s %8s %8s %8d %8d %s (SKIPPED)\n")
        % "N/A"
        % "N/A"
        % "N/A"
        % n
        % traj[i]->nframes()
        % traj[i]->filename();
    else
    {
      cout << boost::format("%5d %8d %8d %8d %8d %s\n")
        % j
        % start_cnt
        % (start_cnt + n - 1)
        % n
        % traj[i]->nframes()
        % traj[i]->filename();
      ++j;
    }
    start_cnt += n;
  }
}


// @endcond



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection* sopts = new opts::BasicSelection("all");
  opts::OutputTrajectoryTypeOptions* otopts = new opts::OutputTrajectoryTypeOptions();
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(otopts).add(topts);
  if (!options.parse(argc, argv)) {
    cerr << "Note- available model file formats (filename suffix) are:\n";
    cerr << availableSystemFileTypes("\t");
    cerr << "Note- available trajectory file formats (filename suffix) are:\n";
    cerr << availableTrajectoryFileTypes("\t");
    exit(-1);
  }


  verbose = bopts->verbosity;
  if (verbose)
    cout << "# " << hdr << endl;

  AtomicGroup model = createSystem(model_name);
  selection = sopts->selection;
  AtomicGroup subset = selectAtoms(model, selection);
  if (subset.empty()) {
    cerr << "Error- no atoms selected in subset\n";
    exit(-10);
  }

  AtomicGroup centered;
  if (!center_selection.empty()) {
    centered = selectAtoms(subset, center_selection);
    if (centered.empty()) {
      cerr << "Error- no atoms selected for centering\n";
      exit(-10);
    }
  }

  AtomicGroup postcentered;
  if (!post_center_selection.empty()) {
    postcentered = selectAtoms(subset, post_center_selection);
    if (postcentered.empty()) {
      cerr << "Error- no atoms selected for post-centering\n";
      exit(-10);
    }
  }

  MultiTrajectory mtraj(traj_names, model, skip, stride);
  if (verbose)
    showTrajectoryTable(mtraj);

  // Wrap since some LOOS tools will expect a pTraj rather than a traj...
  pTraj ptraj(&mtraj, boost::lambda::_1);

  indices = assignTrajectoryFrames(ptraj, topts->range_spec, 0, 1);

  pTrajectoryWriter trajout = otopts->createTrajectory(out_name);
  if (trajout->hasComments())
    trajout->setComments(hdr);

  bool first = true;  // Flag to pick off the first frame for a
                      // reference structure

  // If reimaging, break out the subsets to iterate over...
  vector<AtomicGroup> molecules;
  if (reimage_mode != NONE ) {
    if (!model.hasBonds()) {
      cerr << "WARNING- the model has no connectivity.  Assigning bonds based on distance.\n";
      model.findBonds();
    }

    if (model.hasBonds())
      molecules = model.splitByMolecule();
    else
      molecules = model.splitByUniqueSegid();

    if (verbose)
      cout << boost::format("Reimaging %d molecules\n") % molecules.size();
  }

  // Setup for progress output...
  PercentProgressWithTime watcher;
  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(PercentTrigger(0.25), EstimatingCounter(indices.size()));
  slayer.attach(&watcher);
  if (verbose)
    slayer.start();

  // Iterate over all requested global-frames...
  for (vector<uint>::const_iterator vi = indices.begin(); vi != indices.end(); ++vi) {

    mtraj.readFrame(*vi);
    mtraj.updateGroupCoords(model);

    // Handle Periodic boundary conditions...
    if (box_override) {
      if (first && subset.isPeriodic())
        cerr << "WARNING - overriding existing periodic box.\n";
      model.periodicBox(box);
    }

    // Handle centering...
    if (center_flag) {
      GCoord c = centered.centroid();
      model.translate(-c);
    }


    if (reimage_mode != NONE) {
      if (reimage_mode == AGGRESSIVE || reimage_mode == ZEALOUS) {
        if (reimage_mode == ZEALOUS) {
          for (vGroup::iterator mol = molecules.begin(); mol != molecules.end(); ++mol)
            mol->mergeImage();
        }
        GCoord centroid = centered[0]->coords();
        model.translate(-centroid);
        for (vGroup::iterator mol = molecules.begin(); mol != molecules.end(); ++mol)
          mol->reimage();

        for (uint i=0; i<2; ++i) {
          centroid = centered.centroid();
          model.translate(-centroid);
          for (vGroup::iterator mol = molecules.begin(); mol != molecules.end(); ++mol)
            mol->reimage();
        }

      } else if (reimage_mode == EXTREME) {

        for (vGroup::iterator mol = molecules.begin(); mol != molecules.end(); ++mol) {
          uint midpoint = mol->size() / 2;
          GCoord c = (*mol)[midpoint]->coords();
          mol->translate(-c);
          mol->reimageByAtom();
          mol->translate(c);
        }

        GCoord last_c = centered.centroid();
        bool first = true;
        uint si;
        for (si = 0; si<extreme_max_iters; ++si) {
          GCoord c = centered.centroid();
          if (!first) {
            if (c.distance(last_c) < extreme_threshold)
              break;
          } else
            first = false;
          last_c = c;
          model.translate(-c);
          for (vGroup::iterator mol = molecules.begin(); mol != molecules.end(); ++mol)
            mol->reimage();
        }

        extreme_delta += (last_c.distance(centered.centroid()));
        GCoord c = centered.centroid();
        model.translate(-c);
        extreme_iters += si;

      } else if (reimage_mode == NORMAL){
        for (vGroup::iterator mol = molecules.begin(); mol != molecules.end(); ++mol)
          mol->mergeImage();
      } else {
        cerr << "Error- unknown reimage mode (" << reimage_mode << ") encountered.\n";
        exit(-10);
      }

      if (!post_center_selection.empty()) {
        GCoord postcenter = postcentered.centroid();
        model.translate(-postcenter);
      }

    }

    trajout->writeFrame(subset);

    // Pick off the first frame for the reference structure...
    if (first) {
      PDB pdb = PDB::fromAtomicGroup(subset.copy());
      pdb.remarks().add(hdr);

      if (selection != "all")
        pdb.pruneBonds();

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

  if (reimage_mode == EXTREME && verbose > 2) {
    double avg = static_cast<double>(extreme_iters) / indices.size();
    cerr << boost::format("Average extreme reimage iters = %f\n") % avg;
    cerr << boost::format("Average extreme reimage delta = %f\n") % (extreme_delta / indices.size());
  }
}
