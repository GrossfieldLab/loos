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


using namespace std;
using namespace loos;
namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

typedef vector<AtomicGroup>            vGroup;

string fullHelpMessage(void) {
  string s = 
    "\n"
    " SYNOPSIS\n"
    "\n"
    " Read a single structure and translate a specific selection to the origin\n"
    "\n"
    " DESCRIPTION\n"
    "\n"
    " This tool reads a file with coordinates and translates and reimages it\n"
    " such that a selected set of atoms are at the origin.  The user can\n"
    " specify a subset of atoms to be written out.  The resulting PDB file\n"
    " is written to stdout.\n"
    "\n"
    " Options\n"
    "\n"
    " --center            Selection specifying the atoms that would be moved to the \n"
    "                     origin by the translation.  Defaults to 'all'.\n"
    " --apply             Selection specifying the atoms to which the translation\n"
    "                     is applied.  Defaults to 'all'.\n"
    " --write             Selection specifying which atoms should be be written \n"
    "                     out. Defaults to 'all'.\n"
    " --reimage           If specified, the system will be reimaged, first by \n"
    "                     segment, then by molecule.  If this flag is\n"
    "                     specified, the system must contain connectivity and \n"
    "                     periodicity information.\n"
    " --center_xy         Apply the translation only in the x-y plane\n"
    " --bonds             Specify a second file to use to find the connectivity \n"
    "                     information, e.g. a PSF file to complement a PDB file \n"
    "                     containing cooordinates.\n"
    " --translate=(x,y,z) Translate center to this coordinate\n"
    "\n"
    "\n"
    " EXAMPLES\n"
    "\n"
    "\tcenter-molecule model.pdb >centered.pdb\n"
    "This centers based on all atoms in the model, transforming all atoms\n"
    "\n"
    "\tcenter-molecule --center 'name == \"CA\"' model.pdb >centered.pdb\n"
    "This centers the model based on all alpha-carbons, transforming all atoms.\n"
    "\n"
    "\tcenter-molecule --center 'name == \"CA\"' --write '!hydrogen' model.pdb >centered.pdb\n"
    "This centers the model based on all alpha-carbons, transforming all atoms, but\n"
    "only writing out non-hydrogen atoms.\n"
    "\n"
    "NOTES\n"
    "\n"
    "\tThis tool is largely redundant with other tools, such as \n"
    "\treimage-by-molecule, recenter-trj, and merge-traj.\n"
    "\tUse quotes around the coordinates if you want to translate\n"
    "\n"
    ;
  return(s);
}

// @cond TOOL_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() :
    center_sel("all"),
    apply_sel("all"),
    write_sel("all"),
    bonds_name(""),
    reimage(false),
    center_xy(false)
  { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("center", po::value<string>(&center_sel)->default_value(center_sel), "Selection to calculate the offset from")
      ("apply", po::value<string>(&apply_sel)->default_value(apply_sel), "Selection to actually center")
      ("write", po::value<string>(&write_sel)->default_value(write_sel), "Selection to write to stdout")
      ("reimage", po::value<bool>(&reimage)->default_value(reimage), "Reimage by molecule after")
      ("center_xy", po::value<bool>(&center_xy)->default_value(center_xy), "Center only x&y dimensions")
      ("bonds", po::value<string>(&bonds_name), "Use this model for connectivity")
      ("translate", po::value<GCoord>(&translate)->default_value(GCoord(0.0,0.0,0.0)), "Translate center to this location");
    
  }

  string print() const {
    ostringstream oss;

    oss << boost::format("center='%s',apply='%s',write='%s',reimage=%d,center_xy=%d,bonds='%s',translate=%s")
      % center_sel
      % apply_sel
      % write_sel
      % reimage
      % center_xy
      % bonds_name
      % translate;

    return(oss.str());
  }


  string center_sel, apply_sel, write_sel, bonds_name;
  GCoord translate;
  bool reimage, center_xy;
};
// @endcond




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
  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::ModelWithCoords* mopts = new opts::ModelWithCoords;
  ToolOptions* topts = new ToolOptions;
  opts::AggregateOptions options;
  options.add(bopts).add(mopts).add(topts);

  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup model = mopts->model;
  
  if (topts->reimage) {
    if (!model.isPeriodic()) {
      cerr << "WARNING- Reimaging requested, but the model has no periodic box information\n";
    } else {
      if (!topts->bonds_name.empty()) {
        AtomicGroup bonds = createSystem(topts->bonds_name);
        copyBonds(model, bonds);
      }

      if (!model.hasBonds()) {
        cerr << "WARNING- The model has no connectivity.  Assigning bonds based on distance.\n";
        model.findBonds();
      }
    }
  }

  AtomicGroup center_mol = selectAtoms(model, topts->center_sel);
  GCoord center = center_mol.centroid();
  if (topts->center_xy) center.z() = 0.0;

  AtomicGroup apply_mol = selectAtoms(model, topts->apply_sel);
  GCoord offset = topts->translate - center;
  
  for (AtomicGroup::iterator atom = apply_mol.begin(); atom != apply_mol.end(); ++atom)
    (*atom)->coords() += offset;

  if (topts->reimage && model.isPeriodic()) {
    vGroup molecules = model.splitByMolecule();
    vGroup segments = model.splitByUniqueSegid();
      
    for (vGroup::iterator seg = segments.begin(); seg != segments.end(); ++seg)
      seg->reimage();

    for (vGroup::iterator mol =molecules.begin(); mol != molecules.end(); ++mol)
      mol->reimage();
  }

  AtomicGroup write_mol = selectAtoms(model, topts->write_sel);
  PDB pdb = PDB::fromAtomicGroup(write_mol);
  pdb.remarks().add(hdr);
  cout << pdb;
}
