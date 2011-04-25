/*
  water-sides.cpp
  (c) 2008, 2009 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School


  Desc:
    Classifies a water as being on one side of the membrane or the
    other or inside the membrane (1 = upper, 0 = inside, -1 = lower).

*/


#include <loos.hpp>
#include <limits>
#include <boost/format.hpp>
#include <boost/program_options.hpp>


using namespace std;
using namespace loos;
namespace po = boost::program_options;

typedef std::pair<double,double> Range;
typedef Math::Matrix<int, Math::ColMajor> Matrix;

Range membrane(0.0, 0.0);
string model_name, traj_name, selection_string;

enum Location { UPPER = 1, MEMBRANE = 0, LOWER = -1 };


Range parseRange(const string& s) {
  double a, b;
  int i = sscanf(s.c_str(), "%lf:%lf", &a, &b);
  if (i != 2) {
    cerr << "Parse error with " << s << endl;
    exit(-1);
  }
  return(Range(a,b));
}


void parseOptions(int argc, char *argv[]) {

  try {
    string s;

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("selection,s", po::value<string>(&selection_string)->default_value("name == 'OH2'"), "Atoms to calculate over")
      ("membrane,m", po::value<string>(), "Range for the membrane (zmin:zmax)");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Trajectory filename");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj"))) {
      cerr << "Usage- " << argv[0] << " [options] model-name trajectory-name\n";
      cerr << generic;
      exit(-1);
    }

    if (vm.count("membrane"))
      membrane = parseRange(vm["membrane"].as<string>());
  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  AtomicGroup model = createSystem(model_name);
  pTraj traj = createTrajectory(traj_name, model);
  AtomicGroup subset = selectAtoms(model, selection_string);

  uint m = subset.size();
  uint n = traj->nframes();
  
  Matrix M(m,n);
  for (uint i=0; i<n; ++i) {
    traj->readFrame(i);
    traj->updateGroupCoords(subset);
    for (uint j=0; j<m; ++j) {
      GCoord c = subset[j]->coords();
      Location l;
      if (c[2] > membrane.second)
        l = UPPER;
      else if (c[2] >= membrane.first)
        l = MEMBRANE;
      else
        l = LOWER;
      M(j,i) = l;
    }
  }

  writeAsciiMatrix(cout, M, hdr);
}
