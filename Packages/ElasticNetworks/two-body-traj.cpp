/** @file */

//28July10 - looks at side-nodes.cpp (unreleased loos tool made by Tod)

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010 Tod D. Romo
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


typedef vector<AtomicGroup>   vGroup;

//finds an atom in the psf that matches an atom in the pdb
// and copies it into the AtomicGroup to obtain accurate mass
pAtom findMatch(const pAtom& probe, const AtomicGroup& grp) {
  for (AtomicGroup::const_iterator i = grp.begin(); i != grp.end(); ++i)
    if ((*i)->name() == probe->name() && (*i)->id() == probe->id()
        && (*i)->resname() == probe->resname() && (*i)->resid() == probe->resid()
        && (*i)->segid() == probe->segid())
      return(*i);

  pAtom null;
  return(null);
}


string selection, psf_file;


// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("psf", po::value<string>(&psf_file), "Include a psf file for mass information");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("psf='%s'") % psf_file;
    return(oss.str());
  }
};

// @endcond


int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  ToolOptions* topts = new ToolOptions;
  opts::RequiredArguments* ropts = new opts::RequiredArguments("output", "output-prefix");

  opts::AggregateOptions options;
  // Note slightly unusual order here because we want to keep the
  // output as the first required argument (i.e. appearing before the
  // model & trajectory)
  options.add(bopts).add(sopts).add(ropts).add(tropts).add(topts);

  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup model = tropts->model;
  AtomicGroup subset = selectAtoms(model, sopts->selection);
  pTraj traj = tropts->trajectory;

  string out_name = ropts->value("output");

  DCDWriter dcdout(out_name + ".dcd");
  dcdout.setTitle(hdr);
  bool first = true;


  if (!psf_file.empty()) { //if there is a psf file make a new AtomicGroup
    AtomicGroup structure = createSystem(psf_file);
    //loop over the whole subset, and return a match if i* is found in structure
    //this looks in psf->pAtom for an instance of the specified pdb->pAtom
    //if a match is found then the mass is copied from the psf->pAtom
    //to the subset->mass (i.e. copied into the AtomicGroup) 
    for (AtomicGroup::iterator i = subset.begin(); i != subset.end(); ++i) {
      pAtom match = findMatch(*i, structure); //look in psf for matching atom 
      if (!match) {
        cerr << "ERROR- no match found for atom " << **i << endl; //*(*i)
	//i.e. i is a ptr to a pAtom, which in turn pts to an atom
        exit(-1);
      }
      
      if (!match->checkProperty(Atom::massbit)) {
        cerr << "ERROR- Atom has no mass: " << *match << endl;
        exit(-1);
      }
      //inset mass found in psf file
      (*i)->mass(match->mass());
    }
  }
 
  while (traj->readFrame()){//read next frame, false at end --> so this automatically loops
    traj->updateGroupCoords(model);
    vector<AtomicGroup> residues = subset.splitByResidue();//this makes a separate AG for each residue 
                                                           //that is accessable from by residue[#]
    AtomicGroup cg_sites;
    int currid = model.maxId();
    for (vector<AtomicGroup>::iterator vi = residues.begin(); vi != residues.end(); ++vi) {
      // First, pick off the CA and BB atoms for EACH residue
      AtomicGroup thisResidue = (*vi).select(HeavyAtomSelector());
      AtomicGroup CA = thisResidue.select(AtomNameSelector("CA"));
      AtomicGroup BB = thisResidue.select(BackboneSelector());
      if (CA.empty()) {
	//	cerr << "Error- cannot find CA.\n" << *vi;
	//	exit(-10);
	continue;
      }
      double massholder = 0.0;
      for (int bbi = 0; bbi < BB.size(); ++bbi) {
	massholder += BB[bbi]->mass();      
      }
      //CA should be an AtomicGroup of only one atom, the CA of residue *vi
      CA[0]->occupancy(massholder);
      
      //      CA[0]->occupancy(CA[0]->mass());
      //adds pAtom CA's to AtomicGroup cg_sites
      cg_sites += CA[0];
      
      AtomicGroup sidechain = thisResidue.select(NotSelector(BackboneSelector()));
      if (sidechain.empty()) {
	//	cerr << "Warning- No sidechain atoms for:\n" << *vi;
	continue;
      }
      //Make a new atom, "CGS" and assign it the CMS and 
      //sumed weight of the side chain of residue *vi
      GCoord c = sidechain.centerOfMass();
      pAtom pa(new Atom(++currid, "CGS", c));
      //give pa the same resid, resname, and segid as the current CA 
      //prepare to write it into the pdb with the correct format
      pa->resid(CA[0]->resid());
      pa->resname(CA[0]->resname());
      pa->segid(CA[0]->segid());
      //occupancy column of pdb is used to hold mass 
      //b/c it's higher precission than mass column
      //mass comes from psf if it is entered
      double m = sidechain.totalMass();
      pa->occupancy(m);
      //add pAtom CGS to the AtomicGroup cg_sites it will output in new pdb
      cg_sites += pa;
    }
    AtomicGroup writable = cg_sites.copy();
    writable.renumber();
    
    // Now we have to add the connect records:
    vector<AtomicGroup> fuzzyResidues = writable.splitByResidue();
    //i couldn't figure out how to make this work with ::iterator so i did it my usual way:
    for (uint fe = 0; fe < fuzzyResidues.size(); ++fe) {
      //bond CA to CGS
      if (fuzzyResidues[fe].size() > 1 ){
	fuzzyResidues[fe].getAtom(0)->addBond(fuzzyResidues[fe].getAtom(1));
      }
      //bond CA to next CA
      if(fe < fuzzyResidues.size()-1){
	fuzzyResidues[fe].getAtom(0)->addBond(fuzzyResidues[fe+1].getAtom(0));
      }
    }

    dcdout.writeFrame(writable);
    if(first){ //first structure written to pdb for ref
      PDB pdb = PDB::fromAtomicGroup(writable);
      pdb.remarks().add(hdr);
      string out_pdb_name = out_name + ".pdb";
      ofstream ofs(out_pdb_name.c_str());
      ofs << pdb;
      ofs.close();
      first = false;
    }
  }
}

