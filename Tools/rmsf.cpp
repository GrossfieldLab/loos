/*
  rmsf.cpp

  (c) 2008 Tod D. Romo, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry

  Compute the root mean square fluctuations...

  Usage - rmsf [options] model trajectory
  
  options:
    --selection='selection string'
    --range=start:end

*/


#include <loos.hpp>
#include <cmath>
#include <sstream>

#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include "Statter.hpp"


using namespace std;
using namespace loos;


namespace po = boost::program_options;


string selection, modelname, trajname;
unsigned int start_frame, end_frame;



void parseOptions(int argc, char *argv[]) {

  try {
    
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce this help message")
      ("selection,s", po::value<string>(&selection)->default_value("name == 'CA'"), "Atoms to compute RMSF over")
      ("range,r", po::value<string>(), "Range (start:end) of frames to use")
      ("model", po::value<string>(&modelname))
      ("trajectory", po::value<string>(&trajname));

    po::positional_options_description p;
    p.add("model", 1);
    p.add("trajectory", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("trajectory"))) {
      cerr << "Usage- rmsf [options] model trajectory >output\n";
      cerr << desc;
      exit(-1);
    }
    
    if (vm.count("range")) {
      string rangespec = vm["range"].as<string>();
      int n = sscanf(rangespec.c_str(), "%u:%u", &start_frame, &end_frame);
      if (n != 2) {
        cerr << "Range error - " << rangespec << endl;
        exit(-1);
      }
    }

  }    

  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}




int main(int argc, char *argv[]) {
  
  string hdr = invocationHeader(argc, argv);

  parseOptions(argc, argv);

  cout << "# " << hdr << endl;

  AtomicGroup model = createSystem(modelname);
  pTraj traj = createTrajectory(trajname, model);

  AtomicGroup subset = selectAtoms(model, selection);

  if (start_frame == 0 && end_frame == 0) {
    end_frame = traj->nframes() - 1;
  }

  vector<AtomicGroup> frames;
  for (uint frame = start_frame; frame <= end_frame; frame++) {
    traj->readFrame(frame);
    traj->updateGroupCoords(subset);
    AtomicGroup frame = subset.copy();
    frames.push_back(frame);
  }

  AtomicGroup avg = averageStructure(frames);
  uint n = avg.size();
  uint m = frames.size();

  vector<double> rmsf(n, 0.0);
  for (uint i = 0; i < m; i++)
    for (uint j = 0; j < n; j++) {
      double d = frames[i][j]->coords().distance2(avg[j]->coords());
      rmsf[j] += d;
    }

  for (uint i = 0; i < n; i++)
    rmsf[i] = sqrt(rmsf[i] / m);

  cout << "# atomid\tresid\tRMSF\n";
  for (uint i = 0; i < n; i++)
    cout << boost::format("%10d %6d   %f\n") % avg[i]->id() % avg[i]->resid() % rmsf[i];

}
