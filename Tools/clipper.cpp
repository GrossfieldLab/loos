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
namespace po = boost::program_options;

string model_name;
string selection_name;

vector<GCoord> planes;
bool byresidue = false;



void fullHelp(void) {
  cout << "Sorry, no extra help available at this time\n";
}



void parseOptions(int argc, char *argv[]) {
  vector<string> clips;

  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("fullhelp", "Even more help")
      ("byres,b", "Clip by residue (rather than by atom)")
      ("selection,s", po::value<string>(&selection_name)->default_value("all"), "Selection to apply clipping planes to");


    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("clip", po::value< vector<string> >(&clips), "Clipping planes");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("clip", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || vm.count("fullhelp") || !vm.count("model") || clips.empty() || clips.size() % 3 != 0) {
      cerr << "Usage- " << argv[0] << " [options] model-name (p1) (p2) (p3) [(p1) (p2) (p3) ...]\n";
      cerr << generic;
      if (vm.count("fullhelp"))
        fullHelp();
      exit(-1);
    }

    // Process the clipping plane specifications...
    for (vector<string>::iterator i = clips.begin(); i != clips.end(); ++i) {
      istringstream ss(*i);
      GCoord c;
      if (!(ss >> c)) {
        cerr << "*ERROR* Cannot parse coordinates " << *i << endl;
        exit(-10);
      }
      planes.push_back(c);
    }
    
    if (vm.count("byres"))
      byresidue = true;

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}



int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);
  parseOptions(argc, argv);
  AtomicGroup model = createSystem(model_name);
  AtomicGroup subset = selectAtoms(model, selection_name);

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
        if (byresidue) {
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
  for (AtomicGroup::iterator i = model.begin(); i != model.end(); ++i)
    if (!(*i)->checkProperty(Atom::flagbit))
      clipped.append(*i);

  PDB pdb = PDB::fromAtomicGroup(clipped);
  pdb.remarks().add(hdr);
  cout << pdb;
}
