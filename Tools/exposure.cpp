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

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

typedef vector<AtomicGroup> vGroup;


vector<uint> indices;
double inner_cutoff, outer_cutoff;
string probe_selection;
vector<string> target_selections;
bool symmetry = false;
int verbosity = 0;
bool normalize = true;
bool average = true;



string fullHelpMessage(void) {
  string s = "Examples:\n"
    " * exposure simulation.pdb simulation.dcd 'segid == \"PROT\"'\n"
    "   Computes the solvent exposure for the molecule with segid\n"
    "   \"PROT\".\n"
    "\n"
    " * exposure -P 'segid =~ \"^L\"' simulation.pdb simulation.dcd 'resname == \"HEXO\" && segid == \"P1\"'\n"
    "   Computes the exposure of the residue HEXO with segid P1 to a\n"
    "   lipid membrane (assuming the lipids have segids begining with \"L\".\n"
    "   This could be used to determine the degree of insertion of the\n"
    "   residue into the membrane, for example.\n"
    "\n"
    " * exposure -R1 -P 'segid =~ \"^L\"' simulation.pdb simulation.dcd 'segid == \"P1\"'\n"
    "   Similar to above, except that it averages over the entire peptide\n"
    "   with segid P1 and considers periodic boundaries when determining\n"
    "   which atoms are within the probe shell.\n"
    "\n"
    " * exposure -R1 -I2 -P 'segid != \"BULK\"' simulation.pdb simulation.dcd 'segid == \"P1\"'\n"
    "   Computes the degree to which P1 is buried, i.e. the density of non-\n"
    "   water atoms about P1, excluding any atom that is within 2 A of an atom\n"
    "   in P1.  Also considers periodic boundaries when computing distances.\n"
    "\n"
    " * exposure -P '!(segid == \"BULK\" || segid == \"P1\")' simulation.pdb simulation.dcd 'segid == \"P1\"'\n"
    "   Computes the degree to which P1 is buried, ignoring the atoms from P1.\n"
    "\n"
    " Note: Exposure calculations can be quite lengthy for large systems/trajectories.\n"
    "       you may want to add '&& !hydrogen' to your selections if speed is an issue.\n";

  return(s);
}


// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("normalize,N", po::value<bool>(&normalize)->default_value(true), "Normalize by volume (i.e. output is density)")
      ("average,A", po::value<bool>(&average)->default_value(true), "Average contacts over selection")
      ("probe,P", po::value<string>(&probe_selection)->default_value("segid == 'BULK' && name == 'OH2'"), "Subset to compute exposure against")
      ("inner,I", po::value<double>(&inner_cutoff)->default_value(0.0), "Inner cutoff (ignore atoms closer than this)")
      ("outer,O", po::value<double>(&outer_cutoff)->default_value(5.0), "Outer cutoff (ignore atoms further away than this)")
      ("reimage,R", po::value<bool>(&symmetry)->default_value(true), "Consider symmetry when computing distances");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("normalize=%d, average=%d, probe='%s', inner=%f, outer=%f, reimage=%d")
      % normalize
      % average
      % probe_selection
      % inner_cutoff
      % outer_cutoff
      % symmetry;

    return(oss.str());
  }
};
// @endcond





double density(const AtomicGroup& target, const AtomicGroup& probe, const double inner_radius, const double outer_radius) {
  
  double or2 = outer_radius * outer_radius;
  double ir2 = inner_radius * inner_radius;

  double vol_out = 4.0/3.0 * M_PI * or2 * outer_radius;
  double vol_inn = 4.0/3.0 * M_PI * ir2 * inner_radius;
  double vol = vol_out - vol_inn;

  GCoord box = target.periodicBox();
  ulong contacts = 0;

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
    contacts += probed;
  }

  double dens = static_cast<double>(contacts);
  if (average)
    dens /= target.size();
  if (normalize)
    dens /= vol;
  return(dens);
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  ToolOptions* topts = new ToolOptions;
  opts::RequiredArguments* ropts = new opts::RequiredArguments;
  ropts->addVariableArguments("target", "target-selection");
  
  opts::AggregateOptions options;
  options.add(bopts).add(tropts).add(topts).add(ropts);
  if (!options.parse(argc, argv))
    exit(-1);

  verbosity = bopts->verbosity;
  
  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;
  indices = tropts->frameList();

  target_selections = ropts->variableValues("target");

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

  PercentProgressWithTime watcher;
  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(PercentTrigger(0.1), EstimatingCounter(indices.size()));
  slayer.attach(&watcher);
  if (verbosity)
    slayer.start();

  for (vector<uint>::iterator frame = indices.begin(); frame != indices.end(); ++frame) {
    traj->readFrame(*frame);
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
    if (verbosity)
      slayer.update();
  }

  if (verbosity)
    slayer.finish();
}
