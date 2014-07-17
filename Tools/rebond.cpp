/*
  rebond.cpp

  Rebuild connectivity for a model based on distance
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


using namespace std;
using namespace loos;
namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;



string fullHelpMessage(void) {
  string s =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "\tRebuild connectivity for a model based on distance\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tCreate a pdb with connect records from an arbitrary model.\n"
    "This is done based on distance between 'atoms' in the system,\n"
    "although arbitrary models can be used as input (i.e. coarse-\n"
    "grained models or CA-only models).  This program is intended\n"
    "for visualization, NOT to reconstruct acutal connectivity.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\trebond --full 0 --radius 15 --selection 'name == \"CA\"' model.pdb >network.pdb\n"
    "This is useful for visualizing the ENM connection network.  It finds\n"
    "all connections between all CA atoms within 15 Angstroms of each other.\n"
    "Only the CA atoms and their bonds are output in this case.\n"
    "\n"
    "\trebond --radius 15 --selection 'name == \"CA\"' model.pdb >network.pdb\n"
    "Same as above, but will output the entire model.  Any pre-existing bonds\n"
    "stored in the PDB will be removed and only those bonds between CA atoms\n"
    "will be present.\n"
    "\n"
    "\trebond --radius 4 ca_trace.pdb >model.pdb\n"
    "Given a PDB of only CA atoms, this will connect them back into chains.\n"
    "This is useful with the CA-only PDB output from tools like svd.\n"
    "The radius may need to be tweaked...\n"
    "\n"
    "\trebond --radius 15 --selection 'name = \"CA\" && resid < 10'\\\n"
    "\t  --super 'name == \"CA\"' --tag ENV model.pdb >network.pdb\n"
    "The superset selection and tagging are useful for visualizing the\n"
    "connections between the environment and the subset in a VSA calculation.\n"
    "In this example bonds are only calculated between CAs with resid < 10 \n"
    "and all other CAs.  The atoms that belong to the subset are also tagged\n"
    "with the segid 'ENV'.\n"
    "\n"
    "NOTES\n"
    "\n"
    "\tSome visualization programs, such as VMD, have a hard-coded maximum\n"
    "number of bonds that can be displayed.  This may be lower than the\n"
    "real number of bonds when visualizaing ENM networks.  You will need\n"
    "to either recompile your software, or use one that has larger limits,\n"
    "such as PyMol.\n"
    ;


  return(s);
}



// @cond TOOL_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() :
    append_bonds(false),
    full_model_output(true),
    super("all"),
    segid(""),
    radius(1.25)
  { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("superset", po::value<string>(&super)->default_value(super), "Subset to search for bonds against the selection")
      ("radius", po::value<double>(&radius)->default_value(radius), "Radius cutoff for bonding")
      ("add", po::value<bool>(&append_bonds)->default_value(append_bonds), "Add to existing bonds")
      ("tag", po::value<string>(&segid), "Tag the bound atoms with this segid")
      ("full", po::value<bool>(&full_model_output)->default_value(full_model_output), "Output the entire model (or just the subset if =0)");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("superset='%s', radius=%f, add=%d, tag='%s', full=%d")
      % super % radius % append_bonds % segid % full_model_output;
    return(oss.str());
  }

  bool append_bonds, full_model_output;
  string super, segid;
  double radius;
     
};

// @endcond


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::ModelWithCoords* mopts = new opts::ModelWithCoords;
  ToolOptions* topts = new ToolOptions;
  
  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup model = mopts->model;
  if (!topts->append_bonds)
    model.clearBonds();

  AtomicGroup subset = selectAtoms(model, sopts->selection);
  AtomicGroup superset = selectAtoms(model, topts->super);

  for (AtomicGroup::iterator j = subset.begin(); j != subset.end(); ++j) {
    GCoord c = (*j)->coords();
    if (!topts->segid.empty())
      (*j)->segid(topts->segid);

    for (AtomicGroup::iterator i = superset.begin(); i != superset.end(); ++i) {
         if ((*i)->checkProperty(Atom::indexbit)) {
              if ((*i)->index() == (*j)->index())
                   continue;
         } else {
              if ((*i)->id() == (*j)->id())
                   continue;
         }
              
         double d = c.distance((*i)->coords());
         if (d <= topts->radius) {
              (*j)->addBond(*i);
         }

    }
  }

  if (!topts->segid.empty())
    for (AtomicGroup::iterator i = subset.begin(); i != subset.end(); ++i)
      (*i)->segid(topts->segid);

  PDB pdb;
  if (topts->full_model_output)
    pdb = PDB::fromAtomicGroup(model);
  else
    pdb = PDB::fromAtomicGroup(subset);

  pdb.remarks().add(hdr);
  cout << pdb;
}
