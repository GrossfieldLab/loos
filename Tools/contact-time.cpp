/*
  contact-time

  Computes the number of contacts between a probe group and a set of target groups...
*/


/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, 2010 Tod D. Romo
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
#include <limits>

using namespace std;
using namespace loos;
namespace po = boost::program_options;

typedef vector<AtomicGroup> vGroup;


vector<uint> indices;
double inner_cutoff, outer_cutoff;
string probe_selection;
string model_name, traj_name;
vector<string> target_selections;
bool symmetry;
int verbosity;
bool normalize;
bool max_norm;
bool auto_self;

bool fast_filter;
double fast_pad;



void fullHelp() {
  cout <<
    "* Normalization *\n"
    "Normalization can be performed in two ways: row or column.\n"
    "Row normalization gives the percentage contact between the probe\n"
    "and each target relative to all contacts.  Column normalization\n"
    "gives the percentage contact between the probe and each target\n"
    "relative to the maximum number of contacts against the respective\n"
    "target.\n"
    "\n"
    "* Autoself *\n"
    "The autoself option splits the probe selection into a set of\n"
    "molecules based on segid.  It then computes the contacts between\n"
    "all of these molecules (excluding self-to-self) and includes this\n"
    "as an extra column in the output.  As an example, suppose\n"
    "you have a number of AMLPs in a membrane, each with a different\n"
    "segid (i.e. PE1, PE2, ...) and you want to find the percentage\n"
    "contacts between the AMLPs and PEGL, PGGL, and each other.  The\n"
    "command for this would be:\n"
    "\n"
    "contact-time -a1 model.pdb traj.dcd  'segid =~ \"PE\\d+\"'\\\n"
    "      'resname == \"PEGL\"' and 'resname == \"PGGL\"'\n"
    "\n"
    "This will automatically generate a new set of targets based\n"
    "on the probe selection, splitting them into separate molecules\n"
    "based on their segid.  It then computes the unique pair-wise\n"
    "contacts between each AMLP.  The total number of self-contacts\n"
    "is then included as an extra column in the output.\n"
    "\n"
    "* Fast mode *\n"
    "By default, contact-time uses a distance filter to eliminate\n"
    "target atoms that are too far to be considered when looking\n"
    "at each probe atom.  The padding for the radius used to\n"
    "exclude target atoms can be adjusted with the '-p' option.\n"
    "In the unlikely event the filter causes problems, it can\n"
    "be disabled with '-f0'.\n";
}



void parseOptions(int argc, char *argv[]) {
  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help,h", "Produce this help message")
      ("fullhelp", "Even more help")
      ("verbose,v", po::value<int>(&verbosity)->default_value(1), "Enable verbose output")
      ("rownorm,n", po::value<bool>(&normalize)->default_value(true), "Normalize total # of contacts (across row)")
      ("colnorm,c", po::value<bool>(&max_norm)->default_value(false), "Normalize by max value (down a column)")
      ("inner,i", po::value<double>(&inner_cutoff)->default_value(1.5), "Inner cutoff (ignore atoms closer than this)")
      ("outer,o", po::value<double>(&outer_cutoff)->default_value(2.5), "Outer cutoff (ignore atoms further away than this)")
      ("reimage,R", po::value<bool>(&symmetry)->default_value(true), "Consider symmetry when computing distances")
      ("range,r", po::value< vector<string> >(), "Frames of the DCD to use (in Octave-style ranges)")
      ("autoself,a", po::value<bool>(&auto_self)->default_value(false), "Automatically include self-to-self")
      ("fast,f", po::value<bool>(&fast_filter)->default_value(true), "Use the fast-filter method")
      ("fastpad,p", po::value<double>(&fast_pad)->default_value(1.0), "Padding for the fast-filter method");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Trajectory filename")
      ("probe", po::value<string>(&probe_selection), "Probe selection")
      ("target", po::value< vector<string> >(&target_selections), "Target selections");
    
    po::options_description command_line;
    command_line.add(generic).add(hidden);
    
    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);
    p.add("probe", 1);
    p.add("target", -1);

    
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj") && !target_selections.empty() && vm.count("probe"))) {
      cerr << "Usage- " << argv[0] << " [options] model-name trajectory-name probe target [target ...] >output\n";
      cerr << generic;
      exit(-1);
    }

    if (vm.count("fullhelp")) {
      fullHelp();
      exit(0);
    }

    if (normalize && max_norm) {
      cerr << "Error- cannot use both column and row normalization at the same time\n";
      exit(-1);
    }

    if (vm.count("range")) {
      vector<string> ranges = vm["range"].as< vector<string> >();
      indices = parseRangeList<uint>(ranges);
    }

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}




uint contacts(const AtomicGroup& target, const AtomicGroup& probe, const double inner_radius, const double outer_radius) {
  
  double or2 = outer_radius * outer_radius;
  double ir2 = inner_radius * inner_radius;

  GCoord box = target.periodicBox();
  uint contact = 0;

  for (AtomicGroup::const_iterator j = probe.begin(); j != probe.end(); ++j) {
    GCoord v = (*j)->coords();
    
    for (AtomicGroup::const_iterator i = target.begin(); i != target.end(); ++i) {
      GCoord u = (*i)->coords();
      double d = symmetry ? v.distance2(u, box) : v.distance2(u);
      if (d >= ir2 && d <= or2)
        ++contact;
    }
  }

  return(contact);
}


AtomicGroup pickNearbyAtoms(const AtomicGroup& target, const AtomicGroup& probe, const double radius) {

  GCoord c = probe.centroid();
  GCoord box = probe.periodicBox();
  double maxrad = probe.radius() + radius;
  maxrad *= maxrad;

  
  AtomicGroup nearby;
  nearby.periodicBox(probe.periodicBox());
  for (AtomicGroup::const_iterator i = target.begin(); i != target.end(); ++i) {
    double d = symmetry ? c.distance2((*i)->coords(), box) : c.distance2((*i)->coords());
    if (d <= maxrad)
      nearby.append(*i);
  }

  return(nearby);
}


uint fastContacts(const AtomicGroup& target, const vGroup& probes, const double inner, const double outer) {
  uint total_contacts = 0;

  
  for (vGroup::const_iterator i = probes.begin(); i != probes.end(); ++i) {
    AtomicGroup new_target = pickNearbyAtoms(target, *i, outer+fast_pad);
    uint c = contacts(new_target, *i, inner, outer);
    total_contacts += c;
  }

  return(total_contacts);
}


// Given a vector of groups, compute the number of contacts between
// unique pairs of groups, excluding the self-to-self
//
// Note: this assumes that you want to compare the -ENTIRE- molecule
uint autoSelfContacts(const vGroup& selves, const double inner_radius, const double outer_radius) {
  
  uint total_contact = 0;
  uint n = selves.size();
  for (uint j=0; j<n-1; ++j)
    for (uint i=j+1; i<n; ++i) {
      total_contact += contacts(selves[j], selves[i], inner_radius, outer_radius);
    }

  return(total_contact);
}




// Normalize across a row of contacts against targets.
// Skips the first col [expecting time to be stored there]
// This gives the percentage of all contacts against each target in
// the respective columns...
//
// If the normalization would be 0, then leave the elements unchanged,
// but issue a warning

void rowNormalize(DoubleMatrix& M) {

  for (uint j=0; j<M.rows(); ++j) {
    double sum = 0.0;
    for (uint i=1; i<M.cols(); ++i)
      sum += M(j, i);
    if (sum == 0.0) {
      cerr << "WARNING- zero sum in rowNormalize()\n";
      sum = 1.0;
    }
    for (uint i=1; i<M.cols(); ++i)
      M(j, i) /= sum;
  }
}



// Normalize down a column of contacts by the max contacts.
// Skips the first col [expecting time to be stored there]
// This gives the fractional contacts (vs max) for each target over
// time.
//
// If the normalization would be 0, then leave the elements unchanged,
// but issue a warning

void colNormalize(DoubleMatrix& M) {

  for (uint i=1; i<M.cols(); ++i) {
    double max = -numeric_limits<double>::max();
    for (uint j=0; j<M.rows(); ++j)
      if (M(j, i) > max)
        max = M(j, i);

    if (max == 0.0) {
      cerr << "WARNING- zero max in colNormalize()\n";
      max = 1.0;
    }
    for (uint j=0; j<M.rows(); ++j)
      M(j, i) /= max;
  }
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  AtomicGroup model = createSystem(model_name);
  pTraj traj = createTrajectory(traj_name, model);

  // Handle a request for a subset of frames
  if (indices.empty())
    for (uint i=0; i<traj->nframes(); ++i)
      indices.push_back(i);

  AtomicGroup probe = selectAtoms(model, probe_selection);

  // Build each of the requested targets...
  vGroup targets;
  for (vector<string>::iterator i = target_selections.begin(); i != target_selections.end(); ++i)
    targets.push_back(selectAtoms(model, *i));


  // Size of the output matrix
  uint rows = indices.size();
  uint cols = targets.size() + 1;

  // If comparing self, split apart molecules by unique segids
  vGroup myselves;
  if (auto_self || fast_filter) {
    if (auto_self)
      ++cols;
    myselves = probe.splitByUniqueSegid();
  }

  DoubleMatrix M(rows, cols);
  uint t = 0;

  // Setup our progress counter since this can be a time-consuming
  // program, but only if verbose output is requested.
  PercentProgressWithTime watcher;
  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(PercentTrigger(0.1), EstimatingCounter(indices.size()));
  slayer.attach(&watcher);
  if (verbosity)
    slayer.start();

  for (vector<uint>::iterator frame = indices.begin(); frame != indices.end(); ++frame) {
    traj->readFrame(*frame);
    traj->updateGroupCoords(model);

    if (symmetry && !model.isPeriodic()) {
      cerr << "ERROR - the trajectory must be periodic to use --reimage\n";
      exit(-1);
    }

    M(t, 0) = t;
    for (uint i=0; i<targets.size(); ++i) {
      double d;
      if (fast_filter)
        d = fastContacts(targets[i], myselves, inner_cutoff, outer_cutoff);
      else
        d = contacts(targets[i], probe, inner_cutoff, outer_cutoff);

      M(t, i+1) = d;
    }

    if (auto_self)
      M(t, cols-1) = autoSelfContacts(myselves, inner_cutoff, outer_cutoff);

    ++t;
    if (verbosity)
      slayer.update();
  }

  if (verbosity)
    slayer.finish();

  if (normalize) {
    if (verbosity)
      cerr << "Normalizing across the row...\n";
    rowNormalize(M);

  } else if (max_norm) {
    if (verbosity)
      cerr << "Normalizing by max column value...\n";
    colNormalize(M);

  } else
    if (verbosity)
      cerr << "No normalization.\n";

  writeAsciiMatrix(cout, M, hdr);
}
