/*
  center-molecule
  
  Centers a molecule/system
*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008,2010, Tod D. Romo
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
#include <boost/program_options.hpp>


using namespace std;
using namespace loos;
namespace po = boost::program_options;

typedef vector<AtomicGroup>            vGroup;

string model_name, bonds_name;
string center_sel, apply_sel, write_sel;
bool reimage;
bool center_xy;



void parseOptions(int argc, char *argv[]) {

  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("center,c", po::value<string>(&center_sel)->default_value("all"), "Selection to calculate the offset from")
      ("apply,a", po::value<string>(&apply_sel)->default_value("all"), "Selection to actually center")
      ("write,w", po::value<string>(&write_sel)->default_value("all"), "Selection to write to stdout")
      ("reimage,r", po::value<bool>(&reimage)->default_value(false), "Reimage by molecule after")
      ("center_xy,x", po::value<bool>(&center_xy)->default_value(false), "Center only x&y dimensions")
      ("bonds,b", po::value<string>(&bonds_name), "Use this model for connectivity");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !vm.count("model")) {
      cerr << "Usage- " << argv[0] << " [options] model-name >output.pdb\n";
      cerr << generic;
      exit(-1);
    }

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}


void copyBonds(AtomicGroup& target, AtomicGroup& source) {

  if (target.size() != source.size()) {
    cerr << "ERROR- centering model and connectivity model have different number of atoms.\n";
    exit(-1);
  }

  for (uint i=0; i<target.size(); ++i)
    target[i]->setBonds(source[i]->getBonds());
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  AtomicGroup model = createSystem(model_name);
  
  if (reimage) {
    if (!model.isPeriodic()) {
      cerr << "WARNING- Reimaging requested, but the model has no periodic box information\n";
    } else {
      if (!bonds_name.empty()) {
        AtomicGroup bonds = createSystem(bonds_name);
        copyBonds(model, bonds);
      }

      if (!model.hasBonds()) {
        cerr << "WARNING- The model has no connectivity.  Assigning bonds based on distance.\n";
        model.findBonds();
      }
    }
  }

  AtomicGroup center_mol = selectAtoms(model, center_sel);
  GCoord center = center_mol.centroid();
  if (center_xy) center.z() = 0.0;

  AtomicGroup apply_mol = selectAtoms(model, apply_sel);
  for (AtomicGroup::iterator atom = apply_mol.begin(); atom != apply_mol.end(); ++atom)
    (*atom)->coords() -= center;

  if (reimage) {
    vGroup molecules = model.splitByMolecule();
    vGroup segments = model.splitByUniqueSegid();
      
    for (vGroup::iterator seg = segments.begin(); seg != segments.end(); ++seg)
      seg->reimage();

    for (vGroup::iterator mol =molecules.begin(); mol != molecules.end(); ++mol)
      mol->reimage();
  }

  AtomicGroup write_mol = selectAtoms(model, write_sel);
  PDB pdb = PDB::fromAtomicGroup(write_mol);
  pdb.remarks().add(hdr);
  cout << pdb;
}
