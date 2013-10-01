/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2012, Tod D. Romo
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


// Base class for handling different windowing functions
struct Window {
    Window(const uint window_size) : _window_size(window_size) { }

    virtual double weight(const uint t) const =0;

    // Total weight for window
    double sum() const {
	double s = 0.0;
	for (uint i=0; i<_window_size; ++i)
	    s += weight(i);
	return(s);
    }


    uint _window_size;
};


struct UniformWindow : public Window {
  UniformWindow(const uint n) : Window(n) { }
  
  double weight(const uint t) const { return(1.0); }
};


struct CosineWindow : public Window {
  CosineWindow(const uint n) : Window(n) { }
  double weight(const uint t) const {
    double d = static_cast<double>(t) / _window_size - 0.5;
    return(cos(d * M_PI/2.0));
  }
};


// ----------------------------------------------------------------------------------


void zeroCoords(AtomicGroup& model) {
  for (AtomicGroup::iterator i = model.begin(); i != model.end(); ++i)
    (*i)->coords(GCoord(0,0,0));
}

void addCoords(AtomicGroup& avg, const AtomicGroup model, const double scale) {
  for (uint i=0; i<avg.size(); ++i)
    avg[i]->coords() += scale * model[i]->coords();
}

void divideCoords(AtomicGroup& avg, const double d) {
  for (AtomicGroup::iterator i = avg.begin(); i != avg.end(); ++i)
    (*i)->coords() /= d;
}


// ----------------------------------------------------------------------------------


int main(int argc, char *argv[]) {
  if (argc != 8) {
    cerr << "Usage- simple-smoother output.dcd model traj selection window stride cosine|uniform\n";
    exit(0);
  }

  string hdr = invocationHeader(argc, argv);

  int k = 1;
  string output_name(argv[k++]);
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  AtomicGroup subset = selectAtoms(model, argv[k++]);
  int window_size = strtoul(argv[k++], 0, 10);
  uint stride = strtoul(argv[k++], 0, 10);

  uint starting_frame = window_size/2;
  uint ending_frame = traj->nframes() - window_size;
  uint n = (ending_frame - starting_frame) / stride;

  string kernel(argv[k++]);
  Window* window;
  if (kernel == "cosine")
    window = new CosineWindow(window_size);
  else if (kernel == "uniform")
    window = new UniformWindow(window_size);
  else {
      cerr << "Error- unknown kernel type.  Must be cosine or uniform.\n";
      exit(-1);
  }
  
      

  PDB pdb = PDB::fromAtomicGroup(subset);
  pdb.remarks().add(hdr);
  cout << pdb;


  DCDWriter dcd(output_name);
  dcd.setHeader(subset.size(), n, 1e-3, false);
  dcd.writeHeader();

  AtomicGroup frame = subset.copy();
  double scaling = window->sum();

  for (uint j=starting_frame; j<ending_frame; j += stride) {

    zeroCoords(frame);

    // Now average...
    for (int wi = 0; wi < window_size; ++wi) {
      double scale = window->weight(wi);
      traj->readFrame(j + wi - wi/2);
      traj->updateGroupCoords(subset);
      addCoords(frame, subset, scale);
    }
    divideCoords(frame, scaling);
    dcd.writeFrame(frame);

  }
}
