/*
  rmsd2ref.cpp

  Computes rmsds between a selection and either its average
  conformation or a reference model, optionally aligning the selection.
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
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


// @cond TOOLS_INTERNAL



string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tCompute the RMSD between the simulation and its average or a reference structure\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool writes out a time-series of the RMSD between the simulation and either\n"
    "the average structure or a reference structure.  The average structure is determined\n"
    "by first using an iterative alignment algorithm described in Grossfield, et al.\n"
    "Proteins 67, 31–40 (2007).  When using a reference structure, the trajectory is\n"
    "superimposed using a Kabsch alignment algorithm.  In both cases, the subset of atoms\n"
    "that the RMSD is calculated over can be different from the subset used in the alignment.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\trmsd2ref model.pdb simulation.dcd >rmsd.asc\n"
    "Computes the RMSD to the average structure.  Alpha-carbons are used for the alignment\n"
    "and the RMSD is calculated over all non-hydrogen and non-solvent atoms.\n"
    "\n"
    "\trmsd2ref --align 'name =~ \"^(C|O|N|CA)$\"' --rmsd 'segid == \"PROT\"' \\\n"
    "\t  model.pdb simulation.dcd >rmsd.asc\n"
    "Computes the RMSD to the average structure.  The trajectory is aligned using all backbone\n"
    "atoms and the RMSD is calculated over all atoms in the \"PROT\" segment.\n"
    "\n"
    "\trmsd2ref --align \"`cat active-sel`\" --rmsd \"`cat active-sel`\" \\\n"
    "\t  --target inactive.pdb --talign \"`cat inactive-sel`\" \\\n"
    "\t  --trmsd \"`cat inactive-sel`\" active.pdb active.dcd >rmsd.asc\n"
    "This complex example calculates the RMSD between a simulation of an \"active\" structure\n"
    "and the \"inactive\" model.  The selection for the active atoms is taken from the file\n"
    "'active-sel' and the selection for the inactive atoms is taken from the file 'inactive-sel'.\n"
    "The same set of atoms is used for both aligning and calculating the structures.\n"
    "As an example, 'active-sel' could contain:\n"
    "\t\t(resid >= 35 && resid <= 45) && name == \"CA\"\n"
    "while 'inactive-sel' could contain:\n"
    "\t\t(resid >= 65 && resid <= 75) && name == \"CA\"\n"
    "\n"
    "NOTES\n"
    "\tThe selections used for aligning and calculating RMSD must match both in number of atoms\n"
    "selected and in the sequence of atoms (i.e. the first atom in the --align selection is\n"
    "matched with the first atom in the --talign selection.)\n"
    "\n"
    "SEE ALSO\n"
    "\trmsds\n";

  return(msg);
}


class ToolOptions : public opts::OptionsPackage {
public:
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("align", po::value<string>(&alignment)->default_value("name == 'CA'"), "Align using this selection")
      ("rmsd", po::value<string>(&selection)->default_value("!(hydrogen || segid =~ 'SOLV|BULK')"), "Compute the RMSD over this selection")
      ("target", po::value<string>(&target_name), "Compute RMSD against this reference target (must have coordinates)")
      ("talign", po::value<string>(&target_align)->default_value(""), "Selection for target to use to align (default is to use --align)")
      ("trmsd", po::value<string>(&target_selection)->default_value(""), "Compute the RMSD over this selection for the target (default is to use --rmsd)")
      ("tolerance", po::value<double>(&tol)->default_value(1e-6), "Tolerance to use for iterative alignment");
  }


  string print() const {
    ostringstream oss;
    oss << boost::format("align='%s', iterative=%d, target='%s', talign='%s', trmsd='%s', tolerance=%f, rmsd='%s'")
      % alignment
      % target_name
      % target_align
      % target_selection
      % tol
      % selection;

    return(oss.str());
  }


  bool postConditions(po::variables_map& map) {

    if (!target_name.empty()) {
      if (target_align.empty()) {
        target_align = alignment;
        cerr << "Warning: Using --align selection for target\n";
      }
      
      if (target_selection.empty()) {
        target_selection = selection;
        cerr << "Warning: Using --rmsd selection for target\n";
      }
    }
      
    return(true);
  }

  string alignment;
  string target_name;
  string target_align, target_selection;
  double tol;
  string selection;
};


// @endcond




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  cout << "# " << hdr << endl;

  AtomicGroup molecule = tropts->model;
  pTraj ptraj = tropts->trajectory;
  AtomicGroup subset = selectAtoms(molecule, topts->selection);
  vector<uint> indices = tropts->frameList();

  AtomicGroup target;
  AtomicGroup target_subset;

  if (topts->target_name.empty())
    cerr << boost::format("Computing RMSD vs avg conformation using %d atoms from \"%s\".\n") % subset.size() % topts->selection;
  else {
    target = createSystem(topts->target_name);
    target_subset = selectAtoms(target, topts->target_selection);
    if (target_subset.size() != subset.size()) {
      cerr << boost::format("Error- target selection has %u atoms while trajectory selection has %u.\n") % target_subset.size() % subset.size();
      exit(-1);
    }
    cerr << boost::format("Computing RMSD vs target %s using %d atoms from \"%s\".\n") % topts->target_name % subset.size() % topts->selection;
  }

  vector<AtomicGroup> frames;

  // handle aligning, if requested...
  if (! topts->alignment.empty()) {

    // First, parse the alignment selection and extract the
    // appropriate bits from the trajectory model...
    AtomicGroup align_subset = selectAtoms(molecule, topts->alignment);

    // Iteratively align the trajectory...
    if (topts->target_name.empty()) {
      cerr << boost::format("Aligning using %d atoms from \"%s\".\n") % align_subset.size() % topts->alignment;
      
      boost::tuple<vector<XForm>, greal, int> res = iterativeAlignment(align_subset, ptraj, indices, topts->tol);
      vector<XForm> xforms = boost::get<0>(res);
      
      for (uint i = 0; i < indices.size(); ++i) {
        ptraj->readFrame(indices[i]);
        ptraj->updateGroupCoords(subset);
        subset.applyTransform(xforms[i]);
        AtomicGroup frame = subset.copy();
        frames.push_back(frame);
      }

    } else {   // A target was provided and aligning was requested...
      
      AtomicGroup target_align = selectAtoms(target, topts->target_align);
      cerr << boost::format("Aligning using %d atoms from \"%s\".\n") % target_align.size() % topts->alignment;

      for (uint i=0; i<indices.size(); i++) {
        ptraj->readFrame(indices[i]);
        ptraj->updateGroupCoords(molecule);
        GMatrix M = align_subset.superposition(target_align);
        XForm W(M);
        subset.applyTransform(W);
        AtomicGroup frame = subset.copy();
        frames.push_back(frame);
      }

      
    }
    
  } else {  // No aligning was requested, so simply slurp up the trajectory...

    for (uint i=0; i<indices.size(); i++) {
      ptraj->readFrame(indices[i]);
      ptraj->updateGroupCoords(subset);
      AtomicGroup frame = subset.copy();
      frames.push_back(frame);
    }

  }

  // If no external reference structure was specified, set the target
  // to the average of the trajectory...
  if (topts->target_name.empty()) {
    cerr << "Computing using average structure...\n";
    target = averageStructure(frames);
  } else
    target = target_subset;  // Hack!

  vector<double> rmsds;
  double avg_rmsd = 0.0;

  if (frames[0].size() != target.size()) {
    cerr << "Error - trajectory selection and target selection have differing numbers of atoms.\n";
    exit(-10);
  }

  for (uint i=0; i<frames.size(); i++) {
    double d = target.rmsd(frames[i]);
    rmsds.push_back(d);
    avg_rmsd += d;
  }

  avg_rmsd /= frames.size();
  double std_rmsd = 0.0;
  for (uint i=0; i<rmsds.size(); i++) {
    double d = rmsds[i] - avg_rmsd;
    std_rmsd += d*d;
  }
  std_rmsd /= (rmsds.size() - 1);
  std_rmsd = sqrt(std_rmsd);

  cerr << boost::format("Average RMSD was %.3lf, std RMSD was %.3lf\n") % avg_rmsd % std_rmsd;
  for (uint i=0; i<rmsds.size(); i++)
    cout << i << "\t" << rmsds[i] << endl;

}
