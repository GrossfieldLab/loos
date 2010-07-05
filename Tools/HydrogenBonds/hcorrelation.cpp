/*
  hcorrelation.cpp

  Computes time-correlation for hydrogen bonds...
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



typedef vector<double>     vecDouble;
typedef vector<vecDouble>    vecvecDouble;
typedef vector<string>     vString;

// ---------------  GLOBALS

double length_low, length_high;
double max_angle;
bool use_periodicity;
string donor_selection, acceptor_selection;
string model_name;
vString traj_names;
uint maxrows;

// ---------------






void parseArgs(int argc, char *argv[]) {
  
  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("blow,d", po::value<double>(&length_low)->default_value(1.5), "Low cutoff for bond length")
      ("bhi,D", po::value<double>(&length_high)->default_value(3.0), "High cutoff for bond length")
      ("angle,a", po::value<double>(&max_angle)->default_value(30.0), "Max bond angle deviation from linear")
      ("periodic,p", po::value<bool>(&use_periodicity)->default_value(false), "Use periodic boundary")
      ("clip,c", po::value<uint>(&maxrows)->default_value(0), "Clip size of trajectories to this (0 = auto-size)");


    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<vString>(&traj_names), "Traj filename")
      ("donor", po::value<string>(&donor_selection), "Donor selection")
      ("acceptor", po::value<string>(&acceptor_selection), "Acceptor selection");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("donor", 1);
    p.add("acceptor", 1);
    p.add("model", 1);
    p.add("traj", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
      cout << "Usage- " << argv[0] << " [options] sel-1 sel-2 model traj-1 [traj-2 ...]\n";
      cout << generic;
      exit(0);
    }

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}




vecvecDouble extractCorrelations(const BondMatrix& M) {
  vecvecDouble correlations;

  uint m = M.rows();
  uint n = M.cols();

  for (uint i=0; i<n; ++i) {
    bool is_there_anybody_out_there = false;
    
    for (uint j=0; j<m; ++j)
      if (M(j,i)) {
        is_there_anybody_out_there = true;
        break;
      }

    if (!is_there_anybody_out_there)
      continue;

    vecDouble v;
    for (uint j=0; j<m; ++j)
      v.push_back(M(j, i));

    TimeSeries<double> ts(v);

    TimeSeries<double> tc = ts.correl(m/2);
    vecDouble u;
    for (uint j=0; j<tc.size(); ++j)
      u.push_back(tc[j]);
    correlations.push_back(u);
  }

  return(correlations);
}


vecDouble average(const vecvecDouble& A) {
  uint n = A.size();
  uint m = A[0].size();

  vecDouble avg(m,0.0);
  for (uint j=0; j<m; ++j)
    for (uint i=0; i<n; ++i)
      avg[j] += A[i][j];

  for (uint j=0; j<m; ++j)
    avg[j] /= n;

  return(avg);
}


vecDouble stddev(const vecvecDouble& A, const vecDouble& avg) {
  uint n = A.size();
  uint m = A[0].size();

  vecDouble std(m, 0.0);
  if (n <= 3)
    return(std);

  for (uint j=0; j<m; ++j)
    for (uint i=0; i<n; ++i) {
      double d = A[i][j] - avg[j];
      std[j] += d*d;
    }

  for (uint j=0; j<m; ++j)
    std[j] = sqrt(std[j]/(n-1));

  return(std);
}



uint findMinSize(const AtomicGroup& model, const vString& names) {
  uint n = numeric_limits<uint>::max();
  
  for (vString::const_iterator i = names.begin(); i != names.end(); ++i) {
    pTraj traj = createTrajectory(*i, model);
    if (traj->nframes() < n)
      n = traj->nframes();
  }

  return(n);
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  parseArgs(argc, argv);


  AtomicGroup model = createSystem(model_name);

  SimpleAtom::innerRadius(length_low);
  SimpleAtom::outerRadius(length_high);
  SimpleAtom::maxDeviation(max_angle);
  
  
  SAGroup donors = SimpleAtom::processSelection(donor_selection, model, use_periodicity);
  SAGroup acceptors = SimpleAtom::processSelection(acceptor_selection, model, use_periodicity);


  vecvecDouble correlations;
  back_insert_iterator<vecvecDouble> corr_appender(correlations);


  if (maxrows == 0)
    maxrows = findMinSize(model, traj_names);

  cerr << boost::format("Using %d as row cutoff.\n") % maxrows;


  for (vString::const_iterator ci = traj_names.begin(); ci != traj_names.end(); ++ci) {
    cerr << "Processing " << *ci << endl;
    pTraj traj = createTrajectory(*ci, model);
    
    for (SAGroup::iterator j = donors.begin(); j != donors.end(); ++j) {
      BondMatrix bonds = j->findHydrogenBondsMatrix(acceptors, traj, model);
      
      for (uint i=0; i<bonds.cols(); ++i) {
        bool found = false;
        for (uint j=0; j<bonds.rows(); ++j)
          if (bonds(j, i) != 0) { 
            found = true;
            break;
          }

        if (found) {
          TimeSeries<double> ts;
          for (uint j=0; j<bonds.rows(); ++j)
            ts.push_back(bonds(j, i));
          TimeSeries<double> tcorr = ts.correl(maxrows);
          vecDouble vtmp;
          copy(tcorr.begin(), tcorr.end(), back_inserter(vtmp));
          correlations.push_back(vtmp);
        }

      }

    }

  }



  cerr << boost::format("Found %d time-correlations.\n") % correlations.size();

  vecDouble avg = average(correlations);
  vecDouble std = stddev(correlations, avg);
  for (uint j=0; j<avg.size(); ++j)
    cout << j << "\t" << avg[j] << "\t" << std[j] << endl;
  
}

