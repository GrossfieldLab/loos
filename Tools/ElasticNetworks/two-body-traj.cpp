/** @file */

//28July10 - looks at side-nodes.cpp (unreleased loos tool made by Tod)

#include <loos.hpp>
#include <boost/program_options.hpp>

using namespace std;
using namespace loos;
namespace po = boost::program_options;

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
string model_name, out_name, traj_names;
bool CM_backbone = false; //this should be an additional option


void parseOptions(int argc, char *argv[]) {
  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("selection,s", po::value<string>(&selection)->default_value("all"), "Subset selection")
      ("psf,m", po::value<string>(&psf_file), "Include a psf file for mass information")
      ("CM_backbone,b", po::value<string>(), "Compute center of mass of backbone - default uses alpha carbons");
    
  
    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_names), "Trajectory filename")
      ("out", po::value<string>(&out_name), "Output prefix");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("out", 1);
    p.add("model", 1);
    p.add("traj", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj") && vm.count("out"))) {
      cerr << "Usage- " << argv[0] << " [options] output-prefix model-name trajectory-name\n";
      cerr << generic;
      exit(-1);
    }
    
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
  AtomicGroup subset = selectAtoms(model, selection);
  pTraj traj = createTrajectory(traj_names, model);
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
      AtomicGroup CA = (*vi).select(AtomNameSelector("CA"));
      AtomicGroup BB = (*vi).select(BackboneSelector());
      if (CA.empty()) {
	//	cerr << "Error- cannot find CA.\n" << *vi;
	//	exit(-10);
	continue;
      }
      double massholder = 0.0;
      for (uint bbi = 0; bbi < BB.size(); ++bbi) {
	massholder += BB[bbi]->mass();      
      }
      //CA should be an AtomicGroup of only one atom, the CA of residue *vi
      CA[0]->occupancy(massholder);
      
      //      CA[0]->occupancy(CA[0]->mass());
      //adds pAtom CA's to AtomicGroup cg_sites
      cg_sites += CA[0];
      
      AtomicGroup sidechain = (*vi).select(NotSelector(BackboneSelector()));
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

