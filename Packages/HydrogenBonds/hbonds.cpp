/*
  hbonds

  Find putative hydrogen-bonds based on user-specified criteria (angle and distance)
*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010, Tod D. Romo
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
#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include "hcore.hpp"

using namespace std;
using namespace loos;
namespace po = boost::program_options;


typedef vector<string>    veString;
typedef vector<double>    veDouble;
typedef vector<veDouble>  veveDouble;
typedef vector<uint>      veUint;
typedef vector<veUint>    veveUint;

typedef SimpleAtom        SAtom;
typedef vector<SAtom>     SAGroup;



typedef Math::Matrix<double, Math::RowMajor>    Matrix;


// ---------------  GLOBALS

double length_low, length_high;
double max_angle;
bool use_periodicity;
bool verbose;
bool use_stderr;
vector<string> acceptor_names;
vector<string> acceptor_selections;
string donor_selection;
string model_name;
vector<string> traj_names;

uint skip = 0;


// ---------------


void parseArgs(int argc, char *argv[]) {
  
  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("verbose,v", po::value<bool>(&verbose)->default_value(false), "Verbose output")
      ("stderr,s", po::value<bool>(&use_stderr)->default_value(false), "Report stderr rather than stddev")
      ("blow,d", po::value<double>(&length_low)->default_value(1.5), "Low cutoff for bond length")
      ("bhi,D", po::value<double>(&length_high)->default_value(3.0), "High cutoff for bond length")
      ("angle,a", po::value<double>(&max_angle)->default_value(30.0), "Max bond angle deviation from linear")
      ("periodic,p", po::value<bool>(&use_periodicity)->default_value(false), "Use periodic boundary")
      ("acceptor_name,N", po::value< vector<string> >(&acceptor_names), "Name of an acceptor selection (required)")
      ("acceptor,S", po::value< vector<string> >(&acceptor_selections), "Acceptor selection (required)")
      ("skip,k", po::value<uint>(&skip)->default_value(0), "# of frames to skip from the start of the trajectory");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("donor", po::value<string>(&donor_selection), "Donor selection")
      ("model", po::value<string>(&model_name), "Model filename")
      ("trajs", po::value< vector<string> >(&traj_names), "Traj filename");


    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("donor", 1);
    p.add("model", 1);
    p.add("trajs", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("donor") && vm.count("model") && vm.count("trajs"))) {
      cout << "Usage- " << argv[0] << " [options] donor model traj [traj...]\n";
      cout << generic;
      exit(0);
    }

    if (acceptor_selections.empty()) {
      cerr << "Error- must provide at least one acceptor name and selection.\n";
      exit(-1);
    }
    
    if (acceptor_selections.size() != acceptor_names.size()) {
      cerr << "Error- must provide one name for each acceptor selection.\n";
      exit(-1);
    }


  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}



veDouble rowAverage(const Matrix& M) {
  uint m = M.rows();
  uint n = M.cols();

  veDouble avg;

  for (uint j=0; j<m; ++j) {
    double x = 0.0;
    for (uint i=0; i<n; ++i)
      x += M(j, i);
    avg.push_back(x/n);
  }

  return(avg);
}


veDouble rowStd(const Matrix& M, const veDouble& avg) {
  uint m = M.rows();
  uint n = M.cols();


  if (n < 3)
    return(veDouble(m, 0.0));

  veDouble stddev;

  for (uint j=0; j<m; ++j) {
    double x = 0.0;
    for (uint i=0; i<n; ++i) {
      double d = M(j, i) - avg[j];
      x += d*d;
    }
    stddev.push_back(sqrt(x/(n-1.0)));
  }

  return(stddev);
}



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  parseArgs(argc, argv);


  SimpleAtom::innerRadius(length_low);
  SimpleAtom::outerRadius(length_high);
  SimpleAtom::maxDeviation(max_angle);

  cout << "# " << hdr << endl;

  AtomicGroup model = createSystem(model_name);

  SAGroup donors = SimpleAtom::processSelection(donor_selection, model, use_periodicity);

  vector<SAGroup> acceptors;
  for (uint i=0; i<acceptor_selections.size(); ++i) {
    SAGroup acceptor = SimpleAtom::processSelection(acceptor_selections[i], model, use_periodicity);
    acceptors.push_back(acceptor);
  }
  
  acceptor_names.push_back("Unbound/Other");

  uint n = traj_names.size() * donors.size();
  uint m = acceptor_selections.size();



  Matrix M(m+1, n);

  if (verbose)
    cerr << "Processing- ";

  for (uint k = 0; k<traj_names.size(); ++k) {
    if (verbose)
      cerr << traj_names[k] << " ";

    pTraj traj = createTrajectory(traj_names[k], model);
    if (skip >= traj->nframes()) {
      cerr << boost::format("Error- trajectory '%s' only has %d frames in it, but we are skipping %d frames...\n")
        % traj_names[k]
        % traj->nframes()
        % skip;
      exit(-20);
    }

    BondMatrix B(m, donors.size());

    for (uint t = skip; t<traj->nframes(); ++t) {
      traj->readFrame(t);
      traj->updateGroupCoords(model);

      for (uint i=0; i<donors.size(); ++i) {
        for (uint j=0; j<acceptors.size(); ++j) {
          AtomicGroup found = donors[i].findHydrogenBonds(acceptors[j], true);
          if (! found.empty())
            B(j, i) += 1;
        }
      }
    }

    for (uint i=0; i<donors.size(); ++i) {
      double sum = 0.0;
      for (uint j=0; j<m; ++j) {
        double fraction = static_cast<double>(B(j, i)) / traj->nframes();
        sum += fraction;
        M(j, k * donors.size() + i) = fraction;
      }
      M(m, k * donors.size() + i) = (sum > 1.0) ? 0.0 : (1.0 - sum);
    }

  }

  if (verbose)
    cerr << endl;

  veDouble averages = rowAverage(M);
  veDouble standards = rowStd(M, averages);

  for (uint i=0; i<averages.size(); ++i)
    cout << boost::format("%-3d %-20s %.4f %.4f\n")
      % i
      % acceptor_names[i]
      % averages[i]
      % (standards[i] / (use_stderr ? sqrt(donors.size() * traj_names.size()) : 1.0));

}
