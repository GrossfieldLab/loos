/*
  exposure
  
  (c) 2009 Tod D. Romo, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry

  Computes the degree of exposure of a set of selections over time.
  The exposure can be based on either the density of atoms within a
  sphere radius, or the number of "contacts" between a "probe"
  selection and each target selection...
*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008-2009, Tod D. Romo
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

using namespace std;
using namespace loos;
namespace po = boost::program_options;

typedef vector<AtomicGroup> vGroup;


double inner_cutoff, outer_cutoff;
string probe_selection;
string model_name, traj_name;
vector<string> target_selections;
bool volumetric = false;
bool symmetry = false;




void parseOptions(int argc, char *argv[]) {
  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("probe,p", po::value<string>(&probe_selection)->default_value("segid == 'BULK' && name == 'OH2'"), "Subset to compute exposure against")
      ("inner,i", po::value<double>(&inner_cutoff)->default_value(2.0), "Inner cutoff (ignore atoms closer than this)")
      ("outer,o", po::value<double>(&outer_cutoff)->default_value(5.0), "Outer cutoff (ignore atoms further away than this)")
      ("density,d", "Compute exposure by density rather than contacts")
      ("reimage,r", "Consider symmetry when computing distances");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Trajectory filename")
      ("target", po::value< vector<string> >(&target_selections), "Target selections");
    
    po::options_description command_line;
    command_line.add(generic).add(hidden);
    
    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);
    p.add("target", -1);

    
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj") && !target_selections.empty())) {
      cerr << "Usage- " << argv[0] << " [options] model-name trajectory-name target [target ...]\n";
      cerr << generic;
      exit(-1);
    }

    if (vm.count("density"))
      volumetric = true;
    
    if (vm.count("symmetry"))
      symmetry = true;

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}



double exposure(const AtomicGroup& target, const AtomicGroup& probe, const AtomicGroup& not_probe, const double inner, const double outer) {
  double ic2 = inner * inner;
  double oc2 = outer * outer;
  ulong probed = 0;
  ulong ntotal = 0;
  GCoord box = target.periodicBox();

  for (AtomicGroup::const_iterator j = target.begin(); j != target.end(); ++j) {
    GCoord v = (*j)->coords();

    for (AtomicGroup::const_iterator i = probe.begin(); i != probe.end(); ++i) {
      GCoord u = (*i)->coords();
      double d = symmetry ? v.distance2(u, box) : v.distance2(u);
      if (d >= ic2 && d <= oc2) {
        ++probed;
        ++ntotal;
      }
    }

    for (AtomicGroup::const_iterator i = not_probe.begin(); i != not_probe.end(); ++i) {
      GCoord u = (*i)->coords();
      double d = v.distance2(u);
      if (d >= ic2 && d <= oc2) {
        ++ntotal;
      }
    }
  }

  return(static_cast<double>(probed) / ntotal);
}



double density(const AtomicGroup& target, const AtomicGroup& probe, const double radius) {
  double r2 = radius * radius;

  
  double dens = 0.0;
  double vol = 4.0/3.0 * M_PI * pow(radius, 3.0);
  GCoord box = target.periodicBox();
  
  for (AtomicGroup::const_iterator j = target.begin(); j != target.end(); ++j) {
    GCoord v = (*j)->coords();
    ulong probed = 0;
    
    for (AtomicGroup::const_iterator i = probe.begin(); i != probe.end(); ++i) {
      GCoord u = (*i)->coords();
      double d = symmetry ? v.distance2(u, box) : v.distance2(u);
      if (d <= r2) {
        ++probed;
      }
    }

    double d = static_cast<double>(probed) / vol;
    dens += d;
  }

  dens /= target.size();
  return(dens);
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  AtomicGroup model = createSystem(model_name);
  pTraj traj = createTrajectory(traj_name, model);

  AtomicGroup probe = selectAtoms(model, probe_selection);
  AtomicGroup not_probe = selectAtoms(model, "!(" + probe_selection + ")");

  vGroup targets;
  for (vector<string>::iterator i = target_selections.begin(); i != target_selections.end(); ++i)
    targets.push_back(selectAtoms(model, *i));

  cout << "# " << hdr << endl;
  cout << "# t ";
  for (uint i=0; i < target_selections.size(); ++i)
    cout << (volumetric ? "Density_" : "Ratio_") << i << "\t";
  cout << endl;

  ulong t = 0;

  while (traj->readFrame()) {
    traj->updateGroupCoords(model);

    if (symmetry && !model.isPeriodic()) {
      cerr << "ERROR - the trajectory must be periodic to use --reimage\n";
      exit(-1);
    }

    cout << boost::format("%8d") % t;

    for (vGroup::const_iterator i = targets.begin(); i != targets.end(); ++i) {
      double d;
      if (volumetric)
        d = density(*i, probe, outer_cutoff);
      else
        d = exposure(*i, probe, not_probe, inner_cutoff, outer_cutoff);
      cout << boost::format("  %8.6f") % d;
    }
    cout << endl;
    ++t;
  }
}
