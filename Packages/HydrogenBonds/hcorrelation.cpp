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

using namespace HBonds;
using namespace std;
using namespace loos;
namespace po = boost::program_options;
namespace opts = loos::OptionsFramework;



typedef vector<double>     vecDouble;
typedef vector<vecDouble>    vecvecDouble;
typedef vector<string>     vString;

// ---------------  GLOBALS

double length_low, length_high;
double max_angle;
bool use_periodicity;
bool use_stderr;
string donor_selection, acceptor_selection;
string model_name;
vString traj_names;
uint maxtime;
uint skip;
bool any_hydrogen;

// ---------------




// ---------------
// @cond TOOLS_INTERNAL



string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tHydrogen bond correlation times\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool generates the auto-correlation of putative hydrogen bonds for a trajectory.\n"
    "Given a donor (hydrogen atom) and a set of acceptors, a matrix is constructed for\n"
    "the trajectory where each row corresponds to a frame in the trajectory and each column\n"
    "corresponds to a potential acceptor.  If there is a hydrogen bond present, subject to\n"
    "distance range and angle cutoff, then a 1 is placed in the matrix, otherwise a 0.\n"
    "An auto-correlation is then calculated for each column.  Only columns where there is\n"
    "at least one hydrogen-bond present are included.  If the any-hydrogen flag is set (--any=1),\n"
    "then the state of the hydrogen bond at each time point is the union of all possible\n"
    "acceptors.  This is useful for asking what the correlation is between a donor and -any-\n"
    "possible acceptor.\n"
    "\tThis process is repeated for all possible donors and over all trajectories.  The\n"
    "correlation time-series is then averaged together, so what is written out is the average\n"
    "correlation at a given time, over all donors and all trajectories.  The maximum correlation\n"
    "time is set automatically based on the shortest trajectory.  However, it may be explicitly\n"
    "set with the --maxtime T option.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\thcorrelation 'segid == \"PE1\" && resid == 4 && name == \"HE1\"' \\\n"
    "\t  'name == \"O1\" && (resname == \"PALM\"' model.psf sim.dcd\n"
    "This example uses the HE1 hydrogen of residue 4 in segment PE1 as the donor, and the O1\n"
    "palmitoyl carnonyl oxygen as the acceptor.  The average correlation over all carbonyl\n"
    "oxygens is written out.\n"
    "\n"
    "\thcorrelation --any=1 'segid == \"PE1\" && resid == 4 && name == \"HE1\"'\\\n"
    "\t  'name == \"O1\" && (resname == \"PALM\"' model.psf sim.dcd\n"
    "This example is the same as above, however the correlation is for the peptide hydrogen (HE1)\n"
    "hydrogen bonding to -any- palmitoyl carbonyl.\n"
    "\n"
    "\thcorrelation --blow=2.0 --bhi=4.0 --angle=25.0 --any=1 \\\n"
    "\t  'segid == \"PE1\" && resid == 4 && name == \"HE1\"' \\\n"
    "\t  'name == \"O1\" && (resname == \"PALM\"' model.psf sim.dcd\n"
    "This example is the same as above, but with the hydrogen-bond criteria changed to be\n"
    "2.0 <= distance <= 4.0 and the angle <= 25.0 degrees.\n"
    "\n"
    "SEE ALSO\n"
    "\thmatrix, hbonds\n";

  return(msg);
}


class ToolOptions : public opts::OptionsPackage {
public:
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("blow", po::value<double>(&length_low)->default_value(1.5), "Low cutoff for bond length")
      ("bhi", po::value<double>(&length_high)->default_value(3.0), "High cutoff for bond length")
      ("angle", po::value<double>(&max_angle)->default_value(30.0), "Max bond angle deviation from linear")
      ("periodic", po::value<bool>(&use_periodicity)->default_value(false), "Use periodic boundary")
      ("maxtime", po::value<uint>(&maxtime)->default_value(0), "Max time for correlation (0 = auto-size)")
      ("any", po::value<bool>(&any_hydrogen)->default_value(false), "Correlation for ANY hydrogen bound")
      ("stderr", po::value<bool>(&use_stderr)->default_value(0), "Report standard error rather than standard deviation");

  }

  void addHidden(po::options_description& o) {
    o.add_options()
      ("donor", po::value<string>(&donor_selection), "donor selection")
      ("acceptor", po::value<string>(&acceptor_selection), "acceptor selection")
      ("model", po::value<string>(&model_name), "model")
      ("trajs", po::value< vector<string> >(&traj_names), "Trajectories");
  }

  void addPositional(po::positional_options_description& opts) {
    opts.add("donor", 1);
    opts.add("acceptor", 1);
    opts.add("model", 1);
    opts.add("trajs", -1);
  }

    bool check(po::variables_map& vm) 
        {
            return(
                donor_selection.empty()
                || acceptor_selection.empty()
                || model_name.empty()
                || traj_names.empty()
                   );
        }
    
    
  string help() const {
    return("donor-selection acceptor-selection model traj [traj ...]");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("skip=%d,stderr=%d,blow=%f,bhi=%f,angle=%f,periodic=%d,maxtime=%d,any=%d,acceptor=\"%s\",donor=\"%s\",model=\"%s\",trajs=\"%s\"")
      % skip
      % use_stderr
      % length_low
      % length_high
      % max_angle
      % use_periodicity
      % maxtime
      % any_hydrogen
      % acceptor_selection
      % donor_selection
      % model_name
      % vectorAsStringWithCommas(traj_names);

    return(oss.str());
  }

};




// @endcond




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


  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions opts;
  opts.add(bopts).add(topts);
  if (! opts.parse(argc, argv))
    exit(-1);


  AtomicGroup model = createSystem(model_name);

  SimpleAtom::innerRadius(length_low);
  SimpleAtom::outerRadius(length_high);
  SimpleAtom::maxDeviation(max_angle);
  
  
  SAGroup donors = SimpleAtom::processSelection(donor_selection, model, use_periodicity);
  SAGroup acceptors = SimpleAtom::processSelection(acceptor_selection, model, use_periodicity);


  vecvecDouble correlations;
  back_insert_iterator<vecvecDouble> corr_appender(correlations);


  if (maxtime == 0)
    maxtime = findMinSize(model, traj_names) / 2;

  cerr << boost::format("Using %d as max time for correlation.\n") % maxtime;


  for (vString::const_iterator ci = traj_names.begin(); ci != traj_names.end(); ++ci) {
    cerr << "Processing " << *ci << endl;
    pTraj traj = createTrajectory(*ci, model);
    
    for (SAGroup::iterator j = donors.begin(); j != donors.end(); ++j) {
      if (skip > 0)
        traj->readFrame(skip-1);

      BondMatrix bonds = j->findHydrogenBondsMatrix(acceptors, traj, model);
      if (any_hydrogen) {
        TimeSeries<double> ts;
        for (uint j=0; j<bonds.rows(); ++j) {
          double val = 0.0;
          for (uint i=0; i<bonds.cols(); ++i)
            if (bonds(j, i) != 0) {
              val = 1.0;
              break;
            }
          ts.push_back(val);
        }
        TimeSeries<double> tcorr = ts.correl(maxtime);
        vecDouble vtmp;
        copy(tcorr.begin(), tcorr.end(), back_inserter(vtmp));
        correlations.push_back(vtmp);
            
      } else {
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
            TimeSeries<double> tcorr = ts.correl(maxtime);
            vecDouble vtmp;
            copy(tcorr.begin(), tcorr.end(), back_inserter(vtmp));
            correlations.push_back(vtmp);
          }
          
        }
        
      }

    }

  }



  cerr << boost::format("Found %d time-correlations.\n") % correlations.size();

  vecDouble avg = average(correlations);
  vecDouble std = stddev(correlations, avg);
  double scaling = 1.0;
  if (use_stderr)
    scaling = sqrt(correlations.size());

  cout << "# " << hdr << endl;
  cout << "# Found " << correlations.size() << " time-correlations." << endl;
  cout << "# Using " << ( use_stderr ? "stderr" : "stddev") << endl;

  for (uint j=0; j<avg.size(); ++j)
    cout << j << "\t" << avg[j] << "\t" << (std[j] / scaling) << endl;
  
}

