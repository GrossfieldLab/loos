/*
  rmsds.cpp

  Computes inter-frame rmsds
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

#include <getopt.h>
#include <cstdlib>

typedef unsigned int uint;

struct Globals {
  Globals() : alignment("name == 'CA'"), iterate(false) { }
  
  string alignment;
  bool iterate;
};


Globals globals;


static struct option long_options[] = {
  {"align", required_argument, 0, 'a'},
  {"iterate", no_argument, 0, 'i'},
  {0,0,0,0}
};


static const char* short_options = "a:i";


void show_help(void) {
  Globals defaults;

  cout << "Usage- rmsds [opts] pdb dcd >output\n";
  cout << "       --align=selection    [" << defaults.alignment << "]\n";
  cout << "       --iterate            [" << (defaults.iterate ? string("on") : string("off")) << "]\n";
}


void parseOptions(int argc, char *argv[]) {
  int opt, idx;

  
  while ((opt = getopt_long(argc, argv, short_options, long_options, &idx)) != -1)
    switch(opt) {
    case 'a': globals.alignment = string(optarg); break;
    case 'i': globals.iterate = true; break;
    case 0: break;
    default:
      cerr << "Unknown option '" << (char)opt << "' - ignored.\n";
    }
}


vector<XForm> align(vector<AtomicGroup>& frames, const AtomicGroup& subset, Trajectory& dcd) {
  uint n = dcd.nframes();

  for (uint i = 0; i<n; i++) {
    AtomicGroup frame = subset.copy();
    dcd.readFrame(i);
    dcd.updateGroupCoords(frame);
    frames.push_back(frame);
  }

  boost::tuple<vector<XForm>, greal, int> res = iterativeAlignment(frames, 0.1, 100);
  vector<XForm> xforms = boost::get<0>(res);
  greal rmsd = boost::get<1>(res);
  int iters = boost::get<2>(res);

  cerr << "Subset alignment with " << subset.size()
       << " atoms converged to " << rmsd << " rmsd after "
       << iters << " iterations.\n";

  return(xforms);
}


void readFrames(vector<AtomicGroup>& frames, const AtomicGroup& subset, Trajectory& dcd) {
  uint n = dcd.nframes();

  for (uint i=0; i<n; i++) {
    AtomicGroup frame = subset.copy();
    dcd.readFrame(i);
    dcd.updateGroupCoords(frame);
    frames.push_back(frame);
  }
}


float *interFrameRMSD(vector<AtomicGroup>& frames) {
  uint n = frames.size();
  float *M = new float[n*n];
  uint i;

  for (i=0; i<n*n; i++)
    M[i] = 0.0;

  uint total = n*(n+1)/2;
  uint delta = total / 4;
  uint k = 0;
  uint j;
  for (j=0; j<n; j++) {
    AtomicGroup jframe = frames[j].copy();
    for (i=0; i<=j; i++, k++) {
      double rmsd;

      if (globals.iterate)
	rmsd = jframe.rmsd(frames[i]);
      else {
	(void)jframe.alignOnto(frames[i]);
	rmsd = jframe.rmsd(frames[i]);
      }


      M[j*n+i] = rmsd;
      M[i*n+j] = rmsd;
      
      if (k % delta == 0) {
	float percent = k * 100.0 / total;
	cerr << setprecision(3) << percent << "% complete\n";
      }

    
    }
  }

  return(M);
}


int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  if (argc - optind != 2) {
    cerr << "Invalid arguments.\n";
    show_help();
    exit(-1);
  }

  PDB pdb(argv[optind++]);
  DCD dcd(argv[optind]);

  Parser parsed(globals.alignment);
  KernelSelector selector(parsed.kernel());
  AtomicGroup subset = pdb.select(selector);
  if (subset.size() == 0) {
    cerr << "Error- no atoms selected.\n";
    exit(-1);
  }
  cerr << "Selected " << subset.size() << " atoms.\n";

  vector<AtomicGroup> frames;
  if (globals.iterate) {
    cerr << "Aligning...\n";
    vector<XForm> xforms = align(frames, subset, dcd);
  } else
    readFrames(frames, subset, dcd);

  float *M = interFrameRMSD(frames);


  RawAsciiWriter<float> writer;
  writer.metadata(header);
  writer.write(M, "R", frames.size(), frames.size());
  delete M;
}


  

