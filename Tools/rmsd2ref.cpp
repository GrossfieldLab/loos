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

#include <getopt.h>
#include <cstdlib>

typedef unsigned int uint;

string align_string;
string selection_string;
string model_name, traj_name;
string target_name;

double tolerance = 1e-6;

static struct option long_options[] = {
  {"align", required_argument, 0, 'a'},
  {"target", required_argument, 0, 't'},
  {"tolerance", required_argument, 0, 'T'},
  {"help", no_argument, 0, 'h'},
  {0,0,0,0}
};


static const char* short_options = "a:t:hT:";


void show_help(void) {

  cout << "Computes the RMSD between each frame in a trajectory and a reference structure.\n";
  cout << "  The reference structure can either be an external target or the average conformation.\n";
  cout << "  In each case, you may specify that each frame of the trajectory is aligned against\n";
  cout << "  the target using a different selection.  If no target is specified and aligning is\n";
  cout << "  requested, then an iterative alignment scheme is used to compute a more optimal average\n";
  cout << "  structure.\n\n";
  cout << "Usage- rmsd2ref [options] selection model trajectory\n";
  cout << "       --align=selection_string\n";
  cout << "       --target=model\n";
  cout << "       --tolerance=double\n";
}


void parseOptions(int argc, char *argv[]) {
  int opt, idx;

  
  while ((opt = getopt_long(argc, argv, short_options, long_options, &idx)) != -1)
    switch(opt) {
    case 'a': align_string = string(optarg); break;
    case 't': target_name = string(optarg); break;
    case 'h': show_help(); exit(0);
    case 'T': tolerance = strtod(optarg, 0); break;
    case 0: break;
    default:
      cerr << "Unknown option '" << (char)opt << "' - ignored.\n";
    }

  if (argc - optind != 3) {
    show_help();
    exit(-1);
  }

  selection_string = string(argv[optind++]);
  model_name = string(argv[optind++]);
  traj_name = string(argv[optind++]);
}


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  parseOptions(argc, argv);

  cout << "# " << hdr << endl;

  AtomicGroup molecule = createSystem(model_name);
  pTraj ptraj = createTrajectory(traj_name, molecule);

  cerr << boost::format("Using trajectory \"%s\" with %lu frames.\n") % traj_name % ptraj->nframes();

  Parser parsed(selection_string);
  KernelSelector selector(parsed.kernel());
  AtomicGroup subset = molecule.select(selector);

  AtomicGroup target;
  AtomicGroup target_subset;

  if (target_name.empty())
    cerr << boost::format("Computing RMSD vs avg conformation using %d atoms from \"%s\".\n") % subset.size() % selection_string;
  else {
    target = createSystem(target_name);
    target_subset = target.select(selector);
    if (target_subset.size() != subset.size()) {
      cerr << boost::format("Error- target selection has %u atoms while trajectory selection has %u.\n") % target_subset.size() % subset.size();
      exit(-1);
    }
    cerr << boost::format("Computing RMSD vs target %s using %d atoms from \"%s\".\n") % target_name % subset.size() % selection_string;
  }

  vector<AtomicGroup> frames;

  // handle aligning, if requested...
  if (! align_string.empty()) {

    // First, parse the alignment selection and extract the
    // appropriate bits from the trajectory model...
    Parser parsed_align(align_string);
    KernelSelector align_selector(parsed_align.kernel());
    AtomicGroup align_subset = molecule.select(align_selector);

    // Iteratively align the trajectory...
    if (target_name.empty()) {
      cerr << boost::format("Aligning using %d atoms from \"%s\".\n") % align_subset.size() % align_string;
      
      boost::tuple<vector<XForm>, greal, int> res = loos::iterativeAlignment(align_subset, *ptraj, tolerance);
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
      
      AtomicGroup target_align = target.select(align_selector);
      cerr << boost::format("Aligning using %d atoms from \"%s\".\n") % target_align.size() % align_string;

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
  if (target_name.empty()) {
    cerr << "Computing using average structure...\n";
    target = loos::averageStructure(frames);
  } else
    target = target_subset;  // Hack!

  vector<double> rmsds;
  double avg_rmsd = 0.0;

  assert(frames[0].size() == target.size() && "Error- trajectory selection and target selection have differing number of atoms.");

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
