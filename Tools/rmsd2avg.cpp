/*
  rmsd2avg.cpp

  Computes rmsds between a selection and its average conformation, optionally
  aligning the selection.
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
string pdb_name, dcd_name;

const double tolerance = 1e-6;

static struct option long_options[] = {
  {"align", required_argument, 0, 'a'},
  {0,0,0,0}
};


static const char* short_options = "a:";


void show_help(void) {

  cout << "Usage- rmsd2avg [options] selection pdb dcd\n";
  cout << "       --align=selection_string\n";
}


void parseOptions(int argc, char *argv[]) {
  int opt, idx;

  
  while ((opt = getopt_long(argc, argv, short_options, long_options, &idx)) != -1)
    switch(opt) {
    case 'a': align_string = string(optarg); break;
    case 0: break;
    default:
      cerr << "Unknown option '" << (char)opt << "' - ignored.\n";
    }

  if (argc - optind != 3) {
    show_help();
    exit(-1);
  }

  selection_string = string(argv[optind++]);
  pdb_name = string(argv[optind++]);
  dcd_name = string(argv[optind++]);
}


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  cout << "# " << hdr << endl;

  parseOptions(argc, argv);
  AtomicGroup molecule = createSystem(pdb_name);
  pTraj ptraj = createTrajectory(dcd_name, molecule);

  cerr << boost::format("Using trajectory \"%s\" with %lu frames.\n") % dcd_name % ptraj->nframes();

  Parser parsed(selection_string);
  KernelSelector selector(parsed.kernel());
  AtomicGroup subset = molecule.select(selector);

  cerr << boost::format("Computing RMSD vs avg conformation using %d atoms from \"%s\".\n") % subset.size() % selection_string;

  vector<AtomicGroup> frames;

  if (! align_string.empty()) {
    Parser parsed_align(align_string);
    KernelSelector align_selector(parsed_align.kernel());
    AtomicGroup align_subset = molecule.select(align_selector);

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

  } else {
    uint n = ptraj->nframes();
    for (uint i=0; i<n; i++) {
      ptraj->readFrame(i);
      ptraj->updateGroupCoords(subset);
      AtomicGroup frame = subset.copy();
      frames.push_back(frame);
    }

  }

  AtomicGroup avg = loos::averageStructure(frames);
  vector<double> rmsds;
  double avg_rmsd = 0.0;

  for (uint i=0; i<frames.size(); i++) {
    double d = avg.rmsd(frames[i]);
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
