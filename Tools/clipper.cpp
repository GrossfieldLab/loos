/*
  clipper.cpp

  (c) 2009 Tod D. Romo, Grossfield Lab
  Department of Biochemistry
  University of Rochester School of Medicine and Dentistry

  Applies a set of arbitrary clipping planes to a model, removing
  clipped atoms.
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
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <sstream>

using namespace std;
using namespace loos;
namespace opts = loos::OptionsFramework;


vector<GCoord> planes;

// @cond TOOL_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() : byresidue(false), cliponly(false), auto_selection("") { }

  void addGeneric(opts::po::options_description& o) {
    o.add_options()
      ("byres", opts::po::value<bool>(&byresidue)->default_value(byresidue), "Set to 1 to clip by residue (rather than by atom)")
      ("auto", opts::po::value<string>(&auto_selection)->default_value(auto_selection), "Automatically generate clipping planes for selection")
      ("cliponly", opts::po::value<bool>(&cliponly)->default_value(cliponly), "Set to 1 to only output the clipped selection, not the whole model");
  }

  void addHidden(opts::po::options_description& o) {
    o.add_options()
      ("clip", opts::po::value< vector<string> >(&clips), "Clipping planes");
  }

  void addPositional(opts::po::positional_options_description& pos) {
    pos.add("clip", -1);
  }

  bool check(opts::po::variables_map& map) {
    return( (clips.empty() && auto_selection.empty()) ||
            (clips.size() % 3 != 0) );
  }

  // Clipping planes are specified via 3 points, so must convert them from command-line input
  bool postConditions(opts::po::variables_map& map) {
    for (vector<string>::iterator i = clips.begin(); i != clips.end(); ++i) {
      istringstream ss(*i);
      GCoord c;
      if (!(ss >> c)) {
        cerr << "*ERROR* Cannot parse coordinates " << *i << endl;
        return(false);
      }
      planes.push_back(c);
    }
    return(true);
  }


  bool byresidue, cliponly;
  string auto_selection;
  vector<GCoord> planes;
  vector<string> clips;
};
// @endcond




string fullHelpMessage(void) {
  string msg = 
    "\n"
    "Clipper implements a set of arbitrary clipping planes that can be\n"
    "applied to a selection or to the entire model.  When a selection is\n"
    "used, only the selection is clipped--all other atoms are retained in\n"
    "the output.  Clipping planes are specified by providing three\n"
    "coordinates.  The normal to the plane is determined using the\n"
    "right-hand rule (i.e. assuming the points define the plane in a\n"
    "counter-clockwise fashion).  Atoms that lie on the normal side of the\n"
    "plane are clipped.  Alternatively, if the --byres flag is given, then\n"
    "if an atom is clipped, the entire residue that contains that atom is\n"
    "also clipped regardless of where it lies with respect to the clipping\n"
    "plane.  Finally, any number of clipping planes can be specified on the\n"
    "command line.\n"
    "\n"
    "Examples:\n"
    "\n"
    "  * clipper model.pdb '(0,0,0)' '(1,0,0)' '(0,1,0)'  >clipped.pdb\n"
    "    This defines a clipping plane at z=0 with the normal pointing\n"
    "    along the positive z-axis.\n"
    "\n"
    "  * clipper model.pdb '(0,4,0)' '(1,4,0)' '(0,4,1)'  >clipped.pdb\n"
    "    This defines a clipping plane at y=4 with the normal pointing\n"
    "    along the positive y-axis\n"
    "\n"
    "  * clipper --byres --selection 'segid==\"BULK\"' model.pdb '(0,0,0)' '(1,0,0)' '(0,1,0)'  >clipped.pdb\n"
    "    This defines a clipping plane at z=0 with the normal pointing\n"
    "    along the positive z-axis, but only waters are clipped and if any\n"
    "    water atom is clipped, then the entire water molecule is also\n"
    "    clipped.\n";

  return(msg);
}




void generateClippingPlanes(vector<GCoord>& planes, AtomicGroup& model, const string& sel) {

  AtomicGroup subset = selectAtoms(model, sel);
  vector<GCoord> axes = subset.principalAxes();
  GCoord center = subset.centroid();
  
  planes.push_back(center);
  planes.push_back(center + axes[0]);
  planes.push_back(center + axes[1]);

  cerr << "Automatically adding the following clipping plane:\n\t" << center << endl;
  cerr << "\t" << center + axes[0] << endl;
  cerr << "\t" << center + axes[1] << endl;
}



int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelectionOptions* sopts = new opts::BasicSelectionOptions;
  opts::ModelWithCoordsOptions* mopts = new opts::ModelWithCoordsOptions;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);



  AtomicGroup model = mopts->model;
  AtomicGroup subset = selectAtoms(model, sopts->selection);

  if (!topts->auto_selection.empty())
    generateClippingPlanes(planes, model, topts->auto_selection);
  else
    planes = topts->planes;

  // First, make sure all atoms are unflagged
  for (AtomicGroup::iterator i = model.begin(); i != model.end(); ++i)
    (*i)->clearProperty(Atom::flagbit);

  for (vector<GCoord>::iterator vi = planes.begin(); vi != planes.end();) {
    GCoord x1 = *vi++;
    GCoord x2 = *vi++;
    GCoord x3 = *vi++;
    GCoord n = (x2-x1) ^ (x3-x1);
    n /= n.length();

    for (AtomicGroup::iterator i = subset.begin(); i != subset.end(); ++i) {
      double d = n * ((*i)->coords() - x1);
      if (d >= 0) {
        if (topts->byresidue) {
          AtomicGroup residue = subset.getResidue(*i);
          for (AtomicGroup::iterator j = residue.begin(); j != residue.end(); ++j)
            (*j)->setProperty(Atom::flagbit);
        } else {
          (*i)->setProperty(Atom::flagbit);
        }
      }
    }
  }

  AtomicGroup clipped;
  if (topts->cliponly) {
    for (AtomicGroup::iterator i = subset.begin(); i != subset.end(); ++i)
      if (!(*i)->checkProperty(Atom::flagbit))
        clipped.append(*i);
    
  } else {
    for (AtomicGroup::iterator i = model.begin(); i != model.end(); ++i)
      if (!(*i)->checkProperty(Atom::flagbit))
        clipped.append(*i);
  }

  PDB pdb = PDB::fromAtomicGroup(clipped);
  pdb.remarks().add(hdr);
  cout << pdb;
}
