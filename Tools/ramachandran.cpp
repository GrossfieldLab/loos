/*
  ramachandran.cpp

  
  (c) 2008-2009 Tod D. Romo, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry

  Computes backbone torsion angles for a given set of residues...

  Usage - ramachandran [options] model trajectory selection >output.asc
  
  Notes:

  o The selection must include an extra residue at the beginning and
  end of the segment, i.e. if you want to compute the torsions for
  residues 5-8, then you must use a selection that picks residues 4-9

*/


#include <loos.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <sstream>

using namespace std;
using namespace loos;
namespace po = boost::program_options;


// Convenience typedefs...

typedef vector<AtomicGroup> vGroup;
typedef vector<vGroup> vvGroup;


// ----------------------------------------------------------
// Custom selectors...


// Select an atom based on a specific name
struct AtomNameSelector : public AtomSelector {
  AtomNameSelector(const string& s) : name(s) { }
  bool operator()(const pAtom& pa) const {
    return(pa->name() == name);
  }

  string name;
};





// ----------------------------------------------------------
// Classes to handle extraction of atoms...  This enables easy
// run-time selection of different criteria for which torsions are
// computed... 


// Extractor is the base-class (interface) to how we extract atoms to
// compute torsions for.  It uses NVI to make the extract function
// polymorphic while allowing pre-conditions implemented in the
// base-class.  For more details about NVI,
// see http://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Non-Virtual_Interface
// 


class Extractor {
public:

  // Returns the names of the torsions
  virtual vector<string> names(void) =0;

  // Utility function to make sure each AtomicGroup has only 4 atoms
  void checkFail(const string& s, const AtomicGroup& nascent, const AtomicGroup& residue) {
    if (nascent.size() != 4) {
      cerr << "***ERROR***\n";
      cerr << "Unable to extract atoms for torsion " << s << endl;
      cerr << "Residue dump:\n";
      cerr << residue << endl;
      exit(-10);
    }
  }


  // This is the interface to atom extraction.  Given a vector of
  // AtomicGroups representing residues and an index into the vector,
  // extract the appropriate torsion atoms.

  vGroup extractAtoms(vGroup& residues, int i) {

    // pre-condition making sure the index is valid
    if (i < 1 || i >= static_cast<int>(residues.size())-1)
      throw(runtime_error("Index into residues is out of bounds"));

    return(extractImpl(residues, i));
  }


private:

  // Polymorphic implementation of the actual extraction...  This is
  // what changes if you want to add different types of torsions to
  // use... 
  virtual vGroup extractImpl(vGroup&, int) =0;
};



// Class for extracting phi-psi backbone torsions for proteins...

class PhiPsi : public Extractor {
public:

  // Names, i.e. phi and psi
  vector<string> names(void) {
    vector<string> s;

    s.push_back("phi");
    s.push_back("psi");
    return(s);
  }

  
  // Implementation of function to extract the appropriate atoms for
  // phi & psi...

  vGroup extractImpl(vGroup& residues, int i) {
    
    // Select specific atom types
    AtomNameSelector carbon("C");
    AtomNameSelector nitrogen("N");
    CAlphaSelector calpha;
    
    // Pull out the atoms we'll use.  Since AtomicGroup::select() only
    // returns an AtomicGroup, we use this even though we're actually
    // dealing with a single atom in the group...
    AtomicGroup C_prev = residues[i-1].select(carbon);
    AtomicGroup N = residues[i].select(nitrogen);
    AtomicGroup CA = residues[i].select(calpha);
    AtomicGroup C = residues[i].select(carbon);
    AtomicGroup N_succ = residues[i+1].select(nitrogen);
    
    // Now build an AtomicGroup corresponding to the atoms used to
    // calculate phi
    AtomicGroup phi = C_prev;
    phi.append(N);
    phi.append(CA);
    phi.append(C);
    checkFail("phi", phi, residues[i]);
    
    AtomicGroup psi = N;
    psi.append(CA);
    psi.append(C);
    psi.append(N_succ);
    checkFail("psi", psi, residues[i]);

    // Combine the two AtomicGroup's into a vector & return this...
    vGroup dihe;
    dihe.push_back(phi);
    dihe.push_back(psi);

    return(dihe);
  }
};


// ----------------------------------------------------------

// Class for extracting pseudo-torsions for RNA
// See Wadley, Keating, Duarte, and Pyle (2007) JMB 372:942-957
// for more details...

class PseudoTorsions : public Extractor {
public:
  vector<string> names(void) {
    vector<string> s;

    s.push_back("eta");
    s.push_back("theta");
    return(s);
  }



  vGroup extractImpl(vGroup& residues, int i) {

    // Select specific atom types
    AtomNameSelector C4P("C4'");
    AtomNameSelector P("P");

    // Extract all the atoms we'll use to build the psuedo-torsions
    AtomicGroup c4p_prev = residues[i-1].select(C4P);
    AtomicGroup p = residues[i].select(P);
    AtomicGroup c4p = residues[i].select(C4P);
    AtomicGroup p_succ = residues[i+1].select(P);
    AtomicGroup c4p_succ = residues[i+1].select(C4P);

    // Assemble the group representing eta
    AtomicGroup eta = c4p_prev;
    eta.append(p);
    eta.append(c4p);
    eta.append(p_succ);
    checkFail("eta", eta, residues[i]);

    AtomicGroup theta = p;
    theta.append(c4p);
    theta.append(p_succ);
    theta.append(c4p_succ);
    checkFail("theta", theta, residues[i]);

    // Combine the groups into a vector
    vGroup dihe;
    dihe.push_back(eta);
    dihe.push_back(theta);

    return(dihe);
  }
};

// ----------------------------------------------------------


// Globals...yuck!

vector<uint> indices;    // Frame indices from the trajectory to use

string model_name;       // Name of the model
string traj_name;        // Trajectory
string sel_name;         // Selection

Extractor *extractor;    // Extractor object that's used to specify
                         // which torsions to compute at run-time...


void parseOptions(int argc, char *argv[]) {

  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("pseudo", "Use RNA pseudo-torsions")
      ("range,r", po::value< vector<string> >(), "Frames of the DCD to use (list of Octave-style ranges)");


    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Trajectory filename")
      ("selection", po::value<string>(&sel_name), "Selection");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);
    p.add("selection", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj") && vm.count("selection"))) {
      cerr << "Usage- " << argv[0] << " [options] model-name trajectory-name selection\n";
      cerr << generic;
      exit(-1);
    }

    if (vm.count("range")) {
      vector<string> ranges = vm["range"].as< vector<string> >();
      indices = parseRangeList<uint>(ranges);
    }


    // Instantiate correct Extractor-derived object based on
    // user-selected mode...
    if (vm.count("pseudo"))
      extractor = new PseudoTorsions;
    else
      extractor = new PhiPsi;                // ASSume protein

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}




int main(int argc, char *argv[]) {
  
  string hdr = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  // Read-in/setup objects to access data...
  AtomicGroup model = createSystem(model_name);
  pTraj traj = createTrajectory(traj_name, model);
  AtomicGroup subset = selectAtoms(model, sel_name);

  // If no frame-range specified, use all...
  if (indices.empty()) {
    for (uint i=0; i<traj->nframes(); ++i)
      indices.push_back(i);
  }

  // Split selection into individual residues to operate over...
  vGroup residues = subset.splitByResidue();
  cerr << boost::format("Found %d residues...\n") % residues.size();


  // Data-structure here is a vector of vectors of AtomicGroups.  Each
  // AtomicGroup holds the list of atoms to use for the torsion calc.
  // Each inner vector of AtomicGroups describes the torsions to
  // calculate for each residue/nucleotide.  The outer-vector is the
  // list of all residues/nucleotides to operate over...

  vvGroup torsion_atoms;
  for (int i=1; i<static_cast<int>(residues.size())-1; ++i) {
    vGroup atoms = extractor->extractAtoms(residues, i);
    torsion_atoms.push_back(atoms);
  }
  
  cout << "# " << hdr << endl;
  cout << "# t\t";

  // Construct the header of what torsions were computed...
  vector<string> torsion_names = extractor->names();
  copy(torsion_names.begin(), torsion_names.end(), ostream_iterator<string>(cout, "\t"));
  cout << endl;

  // Iterate over the requested frames from the trajectory...
  vector<uint>::iterator frameno;
  for (frameno = indices.begin(); frameno != indices.end(); ++frameno) {
    traj->readFrame(*frameno);

    // Update ALL atoms...since atoms are shared between AtomicGroups,
    // updating the parent group will update all of the groups we
    // built to calculate torsions...
    traj->updateGroupCoords(model);

    // Iterate over each residue...
    vvGroup::iterator vvi;
    for (vvi = torsion_atoms.begin(); vvi != torsion_atoms.end(); ++vvi) {
      cout << *frameno << " ";

      // Iterate over each group of atoms to use for torsions within
      // the residue...
      vGroup::iterator vi;
      for (vi = (*vvi).begin(); vi != (*vvi).end(); ++vi) {
        double angle = Math::torsion((*vi)[0], (*vi)[1], (*vi)[2], (*vi)[3]);
        cout << angle << "\t";
      }
      cout << endl;
    }
  }

}
