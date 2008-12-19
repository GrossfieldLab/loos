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
#include <boost/format.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

using namespace std;
using namespace loos;



struct Globals {
  string selection;
  string target_name;
  string model_name;
  string traj_name;
  string alignment;
  double tol;
  bool iterate;
};


Globals globals;

void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options", 100);
    generic.add_options()
      ("help", "Produce this help message")
      ("align,a", po::value<string>(&globals.alignment)->default_value("name == 'CA')"), "Align using this selection")
      ("iterative,i", po::value<bool>(&globals.iterate)->default_value(false),"Use iterative alignment method")
      ("target,t", po::value<string>(&globals.target_name), "Compute RMSD against this reference target")
      ("tolerance,T", po::value<double>(&globals.tol)->default_value(1e-6), "Tolerance to use for iterative alignment")
      ("rmsd,r", po::value<string>(&globals.selection)->default_value("!(hydrogen || segid =~ 'SOLV|BULK')"), "Compute the RMSD over this selection");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&globals.model_name), "Model filename")
      ("traj", po::value<string>(&globals.traj_name), "Trajectory filename");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj"))) {
      cerr << "Usage- rmsd2ref [options] model-name trajectory-name\n";
      cerr << generic;
      exit(-1);
    }
  }

  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}
int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  parseOptions(argc, argv);

  cout << "# " << hdr << endl;

  AtomicGroup molecule = createSystem(globals.model_name);
  pTraj ptraj = createTrajectory(globals.traj_name, molecule);

  cerr << boost::format("Using trajectory \"%s\" with %lu frames.\n") % globals.traj_name % ptraj->nframes();

  AtomicGroup subset = selectAtoms(molecule, globals.selection);

  AtomicGroup target;
  AtomicGroup target_subset;

  if (globals.target_name.empty())
    cerr << boost::format("Computing RMSD vs avg conformation using %d atoms from \"%s\".\n") % subset.size() % globals.selection;
  else {
    target = createSystem(globals.target_name);
    target_subset = selectAtoms(target, globals.selection);
    if (target_subset.size() != subset.size()) {
      cerr << boost::format("Error- target selection has %u atoms while trajectory selection has %u.\n") % target_subset.size() % subset.size();
      exit(-1);
    }
    cerr << boost::format("Computing RMSD vs target %s using %d atoms from \"%s\".\n") % globals.target_name % subset.size() % globals.selection;
  }

  vector<AtomicGroup> frames;

  // handle aligning, if requested...
  if (! globals.alignment.empty()) {

    // First, parse the alignment selection and extract the
    // appropriate bits from the trajectory model...
    AtomicGroup align_subset = selectAtoms(molecule, globals.alignment);

    // Iteratively align the trajectory...
    if (globals.target_name.empty()) {
      cerr << boost::format("Aligning using %d atoms from \"%s\".\n") % align_subset.size() % globals.alignment;
      
      boost::tuple<vector<XForm>, greal, int> res = iterativeAlignment(align_subset, ptraj, globals.tol);
      vector<XForm> xforms = boost::get<0>(res);
      
      uint n = ptraj->nframes();
      for (uint i=0; i<n; i++) {
        ptraj->readFrame(i);
        ptraj->updateGroupCoords(subset);
        subset.applyTransform(xforms[i]);
        AtomicGroup frame = subset.copy();
        frames.push_back(frame);
      }

    } else {   // A target was provided and aligning was requested...
      
      AtomicGroup target_align = selectAtoms(target, globals.alignment);
      cerr << boost::format("Aligning using %d atoms from \"%s\".\n") % target_align.size() % globals.alignment;

      uint n = ptraj->nframes();
      for (uint i=0; i<n; i++) {
        ptraj->readFrame(i);
        ptraj->updateGroupCoords(molecule);
        GMatrix M = align_subset.superposition(target_align);
        XForm W(M);
        subset.applyTransform(W);
        AtomicGroup frame = subset.copy();
        frames.push_back(frame);
      }

      
    }
    
  } else {  // No aligning was requested, so simply slurp up the trajectory...

    uint n = ptraj->nframes();
    for (uint i=0; i<n; i++) {
      ptraj->readFrame(i);
      ptraj->updateGroupCoords(subset);
      AtomicGroup frame = subset.copy();
      frames.push_back(frame);
    }

  }

  // If no external reference structure was specified, set the target
  // to the average of the trajectory...
  if (globals.target_name.empty()) {
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
    cout << rmsds[i] << endl;

}
