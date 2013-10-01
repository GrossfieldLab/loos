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


namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

// @cond TOOLS_INTERNAL


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



string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\n"
    "NOTES\n"
      "\n";
  

  return(msg);
}



class ToolOptions : public opts::OptionsPackage {
public:
    ToolOptions() : window_name("cos"),
		    window_size(10),
		    stride(1)
	{ }

    void addGeneric(po::options_description& o) {
        o.add_options()
            ("kernel", po::value<string>(&window_name)->default_value(window_name), "Kernel to use (cos|uniform)")
            ("size", po::value<uint>(&window_size)->default_value(window_size), "Size of window to average over")
            ("stride", po::value<uint>(&stride)->default_value(stride), "How may frames to skip per step");
	
    }


    bool postConditions(po::variables_map& map) 
	{
	    if (window_name == "cos")
		window = new CosineWindow(window_size);
	    else if (window_name == "uniform")
		window = new UniformWindow(window_size);
	    else {
		cerr << "Error- unknown smoothing kernel '" << window_name << "'.\n"
		     << "Must be: cos, uniform\n";
		return(false);
	    }
	    

	    return(true);
	}
    
		
    string print() const {
	ostringstream oss;
	oss << boost::format("kernel='%s',size=%d,stride=%d")
	    % window_name
	    % window_size
	    % stride;
	
	return(oss.str());
    }

    string window_name;
    uint window_size, stride;
    Window* window;
};



// @endcond


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

  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::OutputPrefix* prefopts = new opts::OutputPrefix("smoothed");
  opts::BasicSelection* sopts = new opts::BasicSelection("!hydrogen");
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  ToolOptions* topts = new ToolOptions;
  
  opts::AggregateOptions options;
  options.add(bopts).add(prefopts).add(sopts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
      exit(-1);
  
  string output_name = prefopts->prefix;
  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;
  AtomicGroup subset = selectAtoms(model, sopts->selection);
  int window_size = topts->window_size;
  uint stride = topts->stride;

  uint starting_frame = window_size/2;
  uint ending_frame = traj->nframes() - window_size;
  uint n = (ending_frame - starting_frame) / stride;

  Window* window = topts->window;
      

  PDB pdb = PDB::fromAtomicGroup(subset);
  pdb.remarks().add(hdr);
  string pdb_name = output_name + ".pdb";
  ofstream ofs(pdb_name.c_str());
  ofs << pdb;


  string dcd_name = output_name + ".dcd";
  DCDWriter dcd(dcd_name);
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
