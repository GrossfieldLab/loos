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

using namespace HBonds;
using namespace std;
using namespace loos;
namespace po = boost::program_options;
namespace opts = loos::OptionsFramework;


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
// @cond TOOLS_INTERNAL


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tHydrogen bond occupancy for a trajectory\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool computes the occupancy of putative hydrogen bonds (defined by\n"
    "a simple distance and angle criteria).  The 'donor' selection must have one\n"
    "hydrogen present and the 'acceptor' should have no hydrogens.  Multiple acceptors\n"
    "may be given on the command line.  These are specified by using multiple sets of\n"
    "options, i.e. -N name -S selection where name is the label for the acceptor and\n"
    "selection is the corresponding LOOS selection string.  There must be at least\n"
    "one name/selection pair.  The occupancy calculation can also be performed over\n"
    "multiple trajectories by specifying more than one on the command line.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\thbonds -N 'Carbonyl' -S 'name == \"O1\" && resname == \"PALM\"' \\\n"
    "\t  'resid == 4 && name == \"HE1\"' model.psf traj.dcd\n"
    "This example uses the palmitoyl carbonyl oxygen as the acceptor and the HE1 hydrogen from\n"
    "residue 4 as the donor.\n"
    "\n"
    "\thbonds -N 'Carbonyl' -S 'name == \"O1\" && resname == \"PALM\"' \\\n"
    "\t  -N 'Phosphate' -S 'name == \"OP1\" && resname == \"PEGL\"' \\\n"
    "\t  'resid == 4 && name == \"HE1\"' model.psf traj.dcd\n"
    "This example uses the palmitoyl carbonyl oxygen as above, but also looks for hydrogen\n"
    "bonds with the OP1 phosphate oxygen in residue PEGL.  The same donor as above is used.\n"
    "\n"
    "\thbonds --blow 2 --bhi 4 --angle 20 -N 'Carbonyl' \\\n"
    "\t  -S 'name == \"O1\" && resname == \"PALM\"' 'resid == 4 && name == \"HE1\"' \\\n"
    "\t  model.psf traj.dcd\n"
    "This example is the same as the first, however the criteria for hydrogen bonds are now\n"
    "that they cannot be shorter than 2 angstroms nor longer than 4 angstroms, and the angle\n"
    "cannot be more than 20 degrees from linear.\n"
    "\n"
    "SEE ALSO\n"
    "\thmatrix, hcorrelation\n";

  return(msg);
}



class ToolOptions : public opts::OptionsPackage {
public:
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("skip,k", po::value<uint>(&skip)->default_value(0), "Number of frames to skip")
      ("stderr", po::value<bool>(&use_stderr)->default_value(false), "Report stderr rather than stddev")
      ("blow", po::value<double>(&length_low)->default_value(1.5), "Low cutoff for bond length")
      ("bhi", po::value<double>(&length_high)->default_value(3.0), "High cutoff for bond length")
      ("angle", po::value<double>(&max_angle)->default_value(30.0), "Max bond angle deviation from linear")
      ("periodic", po::value<bool>(&use_periodicity)->default_value(false), "Use periodic boundary")
      ("name,N", po::value< vector<string> >(&acceptor_names), "Name of an acceptor selection (required)")
      ("acceptor,S", po::value< vector<string> >(&acceptor_selections), "Acceptor selection (required)");
  }

  void addHidden(po::options_description& o) {
    o.add_options()
      ("donor", po::value<string>(&donor_selection), "donor selection")
      ("model", po::value<string>(&model_name), "model")
      ("trajs", po::value< vector<string> >(&traj_names), "Trajectories");
  }

  void addPositional(po::positional_options_description& opts) {
    opts.add("donor", 1);
    opts.add("model", 1);
    opts.add("trajs", -1);
  }

  bool postConditions(po::variables_map& map) {
    if (acceptor_selections.empty()) {
      cerr << "Error- must provide at least one acceptor name and selection.\n";
      return(false);
    }
    if (acceptor_selections.size() != acceptor_names.size()) {
      cerr << "Error- must provide one name for each acceptor selection.\n";
      return(false);
    }
    
    return(true);
  }

  string help() const {
    return("donor model traj [traj ...]");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("skip=%d,stderr=%d,blow=%f,bhi=%f,angle=%f,periodic=%d,names=\"%s\",acceptors=\"%s\",donor=\"%s\",model=\"%s\",trajs=\"%s\"")
      % skip
      % use_stderr
      % length_low
      % length_high
      % max_angle
      % use_periodicity
      % vectorAsStringWithCommas(acceptor_names)
      % vectorAsStringWithCommas(acceptor_selections)
      % donor_selection
      % model_name
      % vectorAsStringWithCommas(traj_names);

    return(oss.str());
  }

};




// @endcond


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

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  SimpleAtom::innerRadius(length_low);
  SimpleAtom::outerRadius(length_high);
  SimpleAtom::maxDeviation(max_angle);

  cout << "# " << hdr << endl;
  cout << "# " << vectorAsStringWithCommas(options.print()) << endl;

  AtomicGroup model = createSystem(model_name);

  SAGroup donors = SimpleAtom::processSelection(donor_selection, model, use_periodicity);

  vector<SAGroup> acceptors;
  for (uint i=0; i<acceptor_selections.size(); ++i) {
    SAGroup acceptor = SimpleAtom::processSelection(acceptor_selections[i], model, use_periodicity);
    cout << boost::format("# Group %d size is %d\n") % i % acceptor.size();
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
