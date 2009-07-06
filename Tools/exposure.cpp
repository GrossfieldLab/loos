/*
  exposure
  
  (c) 2009 Tod D. Romo, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry

  Computes the degree of exposure of a set of selections over time.
  The exposure is defined as the average density of a probe selection
  within a spherical shell about each target atom.
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
bool symmetry = false;



void fullHelp(void) {
  cout << "Examples:\n"
    " * exposure simulation.pdb simulation.dcd 'segid == \"PROT\"'\n"
    "   Computes the solvent exposure for the molecule with segid\n"
    "   \"PROT\".\n"
    "\n"
    " * exposure -p 'segid =~ \"^L\"' simulation.pdb simulation.dcd 'resname == \"HEXO\" && segid == \"P1\"'\n"
    "   Computes the exposure of the residue HEXO with segid P1 to a\n"
    "   lipid membrane (assuming the lipids have segids begining with \"L\".\n"
    "   This could be used to determine the degree of insertion of the\n"
    "   residue into the membrane, for example.\n"
    "\n"
    " * exposure -R -p 'segid =~ \"^L\"' simulation.pdb simulation.dcd 'segid == \"P1\"'\n"
    "   Similar to above, except that it averages over the entire peptide\n"
    "   with segid P1 and considers periodic boundaries when determining\n"
    "   which atoms are within the probe shell.\n"
    "\n"
    " * exposure -R -i 2 -p 'segid != \"BULK\"' simulation.pdb simulation.dcd 'segid == \"P1\"'\n"
    "   Computes the degree to which P1 is buried, i.e. the density of non-\n"
    "   water atoms about P1, excluding any atom that is within 2 A of an atom\n"
    "   in P1.  Also considers periodic boundaries when computing distances.\n";
}



void parseOptions(int argc, char *argv[]) {
  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("fullhelp", "Even more help")
      ("probe,p", po::value<string>(&probe_selection)->default_value("segid == 'BULK' && name == 'OH2'"), "Subset to compute exposure against")
      ("inner,i", po::value<double>(&inner_cutoff)->default_value(0.0), "Inner cutoff (ignore atoms closer than this)")
      ("outer,o", po::value<double>(&outer_cutoff)->default_value(5.0), "Outer cutoff (ignore atoms further away than this)")
      ("reimage,R", "Consider symmetry when computing distances");

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

    if (vm.count("help") || vm.count("fullhelp") || !(vm.count("model") && vm.count("traj") && !target_selections.empty())) {
      cerr << "Usage- " << argv[0] << " [options] model-name trajectory-name target [target ...]\n";
      cerr << generic;
      if (vm.count("fullhelp"))
        fullHelp();
      exit(-1);
    }

    if (vm.count("reimage"))
      symmetry = true;

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}




double density(const AtomicGroup& target, const AtomicGroup& probe, const double inner_radius, const double outer_radius) {
  
  double dens = 0.0;
  double or2 = outer_radius * outer_radius;
  double ir2 = inner_radius * inner_radius;

  double vol_out = 4.0/3.0 * M_PI * or2 * outer_radius;
  double vol_inn = 4.0/3.0 * M_PI * ir2 * inner_radius;
  double vol = vol_out - vol_inn;

  GCoord box = target.periodicBox();
  
  for (AtomicGroup::const_iterator j = target.begin(); j != target.end(); ++j) {
    GCoord v = (*j)->coords();
    ulong probed = 0;
    
    for (AtomicGroup::const_iterator i = probe.begin(); i != probe.end(); ++i) {
      GCoord u = (*i)->coords();
      double d = symmetry ? v.distance2(u, box) : v.distance2(u);
      if (d >= ir2 && d <= or2) {
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

  vGroup targets;
  for (vector<string>::iterator i = target_selections.begin(); i != target_selections.end(); ++i)
    targets.push_back(selectAtoms(model, *i));

  cout << "# " << hdr << endl;
  cout << "# t ";
  for (uint i=0; i < target_selections.size(); ++i)
    cout << "Density_" << i << "\t";
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
      d = density(*i, probe, inner_cutoff, outer_cutoff);
      cout << boost::format("  %8.6f") % d;
    }
    cout << endl;
    ++t;
  }
}
