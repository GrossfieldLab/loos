/*
  hmatrix.cpp

  Writes out a matrix representing hbonds over time
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
namespace opts = loos::OptionsFramework;


// ---------------  GLOBALS

double length_low, length_high;
double max_angle;
bool use_periodicity;
string donor_selection, acceptor_selection;
string model_name;
string traj_name;

uint currentTimeStep = 0;


// ---------------


// ---------------
// @cond TOOLS_INTERNAL


class ToolOptions : public opts::OptionsPackage {
public:
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("blow", po::value<double>(&length_low)->default_value(1.5), "Low cutoff for bond length")
      ("bhi", po::value<double>(&length_high)->default_value(3.0), "High cutoff for bond length")
      ("angle", po::value<double>(&max_angle)->default_value(30.0), "Max bond angle deviation from linear")
      ("periodic", po::value<bool>(&use_periodicity)->default_value(false), "Use periodic boundary");
  }

  void addHidden(po::options_description& o) {
    o.add_options()
      ("donor", po::value<string>(&donor_selection), "donor selection")
      ("acceptor", po::value<string>(&acceptor_selection), "acceptor selection");

  }

  void addPositional(po::positional_options_description& opts) {
    opts.add("donor", 1);
    opts.add("acceptor", 1);
  }


  string help() const {
    return("donor-selection acceptor-selection");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("blow=%f,bhi=%f,angle=%f,periodic=%d,acceptor=\"%s\",donor=\"%s\"")
      % length_low
      % length_high
      % max_angle
      % use_periodicity
      % acceptor_selection
      % donor_selection;

    return(oss.str());
  }

};




// @endcond





int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);


  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  ToolOptions* topts = new ToolOptions;
  
  opts::AggregateOptions options;
  options.add(bopts).add(tropts).add(topts);
  if (! options.parse(argc, argv))
    exit(-1);


  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;

  if (use_periodicity && !traj->hasPeriodicBox()) {
    cerr << "Error- trajectory has no periodic box information\n";
    exit(-1);
  }

  SimpleAtom::innerRadius(length_low);
  SimpleAtom::outerRadius(length_high);
  SimpleAtom::maxDeviation(max_angle);

  SAGroup donors = SimpleAtom::processSelection(donor_selection, model, use_periodicity);
  if (donors.size() != 1) {
    cerr << "Error- only specify one donor atom (the attached hydrogen)\n";
    exit(-1);
  }

  SAGroup acceptors = SimpleAtom::processSelection(acceptor_selection, model, use_periodicity);
  BondMatrix bonds = donors[0].findHydrogenBondsMatrix(acceptors, traj, model);
  writeAsciiMatrix(cout, bonds, hdr);
}

