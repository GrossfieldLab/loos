/*
  contact-time
  
  (c) 2010 Tod D. Romo, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry

  Computes the number of contacts between a probe group and a set of target groups...
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
bool symmetry = false;
int verbosity = 0;
bool normalize = true;
bool max_norm = false;
bool local_normal = false;
bool auto_self = false;



void parseOptions(int argc, char *argv[]) {
  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help,h", "Produce this help message")
      ("fullhelp", "Even more help")
      ("verbose,v", "Verbose output")
      ("normalize,n", po::value<bool>(&normalize)->default_value(true), "Normalize total # of contacts")
      ("max,m", po::value<bool>(&max_norm)->default_value(false), "Normalize by max value down a column")
      ("local,l", po::value<bool>(&local_normal)->default_value(false), "Normalize by possible # of contacts (i.e. size of probe selection)")
      ("inner,i", po::value<double>(&inner_cutoff)->default_value(1.5), "Inner cutoff (ignore atoms closer than this)")
      ("outer,o", po::value<double>(&outer_cutoff)->default_value(2.5), "Outer cutoff (ignore atoms further away than this)")
      ("reimage,R", po::value<bool>(&symmetry)->default_value(true), "Consider symmetry when computing distances")
      ("range,r", po::value< vector<string> >(), "Frames of the DCD to use (in Octave-style ranges)")
      ("autoself,a", po::value<bool>(&auto_self)->default_value(false), "Automatically include self-to-self");

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
      cerr << "Usage- " << argv[0] << " [options] model-name trajectory-name probe target [target ...]\n";
      cerr << generic;
      exit(-1);
    }

    if (vm.count("verbose"))
      verbosity = 1;

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




double contacts(const AtomicGroup& target, const AtomicGroup& probe, const double inner_radius, const double outer_radius) {
  
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

  double val = static_cast<double>(contact);
  if (local_normal)
    val /= probe.size();

  return(val);
}



// This assumes that you want to compare the -ENTIRE- molecule
double autoSelfContacts(const vGroup& memyselfandi, const double inner_radius, const double outer_radius) {
  
  double total_contact = 0.0;
  uint n = memyselfandi.size();
  for (uint j=0; j<n-1; ++j)
    for (uint i=j+1; i<n; ++i) {
      total_contact += contacts(memyselfandi[j], memyselfandi[i], inner_radius, outer_radius);
    }

  // Dubious usage...
  if (local_normal)
    total_contact /= memyselfandi[0].size();

  return(total_contact);
}



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

  if (indices.empty())
    for (uint i=0; i<traj->nframes(); ++i)
      indices.push_back(i);

  AtomicGroup probe = selectAtoms(model, probe_selection);

  vGroup targets;
  for (vector<string>::iterator i = target_selections.begin(); i != target_selections.end(); ++i)
    targets.push_back(selectAtoms(model, *i));

  vGroup myself;
  uint rows = indices.size();
  uint cols = targets.size() + 1;
  if (auto_self) {
    ++cols;
    myself = probe.splitByUniqueSegid();
  }
    

  DoubleMatrix M(rows, cols);

  uint t = 0;

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
      d = contacts(targets[i], probe, inner_cutoff, outer_cutoff);
      M(t, i+1) = d;
    }

    if (auto_self)
      M(t, cols-1) = autoSelfContacts(myself, inner_cutoff, outer_cutoff);

    ++t;
    if (verbosity)
      slayer.update();
  }

  if (verbosity)
    slayer.finish();

  if (normalize) {
    cerr << "Normalizing across the row...\n";
    rowNormalize(M);
  } else if (max_norm) {
    cerr << "Normalizing by max column value...\n";
    colNormalize(M);
  } else
    cerr << "No normalization.\n";

  writeAsciiMatrix(cout, M, hdr);
}
