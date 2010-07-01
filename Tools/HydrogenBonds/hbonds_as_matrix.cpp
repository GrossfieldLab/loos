/*
  hbonds_as_matrix.cpp
  (c) 2010 Tod D. Romo, Grossfield Lab, URMC
  

*/


#include <loos.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include "hbonds_core.hpp"

using namespace std;
using namespace loos;
namespace po = boost::program_options;


// ---------------  GLOBALS

double length_low, length_high;
double max_angle;
bool use_periodicity;
string donor_selection, acceptor_selection;
string model_name;
string traj_name;

uint currentTimeStep = 0;


// ---------------



void parseArgs(int argc, char *argv[]) {
  
  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("blow,d", po::value<double>(&length_low)->default_value(1.5), "Low cutoff for bond length")
      ("bhi,D", po::value<double>(&length_high)->default_value(3.0), "High cutoff for bond length")
      ("angle,a", po::value<double>(&max_angle)->default_value(30.0), "Max bond angle deviation from linear")
      ("periodic,p", po::value<bool>(&use_periodicity)->default_value(false), "Use periodic boundary");


    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Traj filename")
      ("donor", po::value<string>(&donor_selection), "Donor selection")
      ("acceptor", po::value<string>(&acceptor_selection), "Acceptor selection");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);
    p.add("donor", 1);
    p.add("acceptor", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
      cout << "Usage- " << argv[0] << " [options] model traj sel-1 sel-2\n";
      cout << generic;
      exit(0);
    }

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  parseArgs(argc, argv);


  AtomicGroup model = createSystem(model_name);
  pTraj traj = createTrajectory(traj_name, model);
  if (use_periodicity && !traj->hasPeriodicBox()) {
    cerr << "Error- trajectory has no periodic box information\n";
    exit(-1);
  }


  SimpleAtom::innerRadius(length_low);
  SimpleAtom::outerRadius(length_high);
  SimpleAtom::maxDeviation(max_angle);

  cout << "# " << hdr << endl;


  SAGroup donors = SimpleAtom::processSelection(donor_selection, model, use_periodicity);
  if (donors.size() != 1) {
    cerr << "Error- only specify one donor atom (the attached hydrogen)\n";
    exit(-1);
  }

  SAGroup acceptors = SimpleAtom::processSelection(acceptor_selection, model, use_periodicity);
  BondMatrix bonds = donors[0].findHydrogenBondsMatrix(acceptors, traj, model);
  writeAsciiMatrix(cout, bonds, hdr);
}

