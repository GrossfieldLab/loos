/*
  ramachandran.cpp

  
  (c) 2008-2009 Tod D. Romo, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry

  Computes backbone torsion angles for a given set of residues...

  Usage - ramachandran [options] model trajectory selection >output.asc
  
  Notes:

  o Not all torsions for a selection can be computed, such as phi-psi
    at the ends of a segment.

  o Missing torsions are replaced with a special value (default is
    -9999)

  o Use the "--skip" flag to exclude residues for which not all
    torsions can be calculated...note that this requires you to select
    one extra residue at either end of the segment

  o Use the "--warn" flag to have ramachandran write out debugging
    info for any residue for which it cannot compute all torsions

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


// Convenience typedefs...

typedef vector<AtomicGroup> vGroup;
typedef vector<vGroup> vvGroup;



// ----------------------------------------------------------
// Classes to handle extraction of atoms...  This enables easy
// run-time selection of different criteria for which torsions are
// computed... 


// Extractor is the base-class (interface) for how we extract atoms to
// compute torsions for.  It uses NVI to make the extract function
// polymorphic while allowing pre- and post-conditions implemented in
// the base-class.  For more details about NVI,
// see http://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Non-Virtual_Interface
// 


class Extractor {
public:

  Extractor() : missing_atoms_warn(false), skip_when_missing(false) { }
  virtual ~Extractor() { }



  // Configuration...
  // Warn if atoms are missing and a torsion cannot be calculated...
  void missingAtomsWarn() { missing_atoms_warn = true; }

  // If atoms are missing, return an empty vGroup so this residue will
  // be skipped in future calculations...
  void skipMissingResidues() { skip_when_missing = true; }

  // Returns the names of the torsions
  virtual vector<string> names(void) =0;

  // Utility function to make sure each AtomicGroup has only 4 atoms
  // for a torsion
  void checkAtoms(const string& s, const AtomicGroup& nascent, const AtomicGroup& residue) {
    if (!missing_atoms_warn)
      return;

    if (nascent.size() != 4) {
      cerr << "***WARNING***\n";
      cerr << "Unable to extract atoms for torsion " << s << endl;
      cerr << "Residue dump:\n";
      cerr << residue << endl;
      cerr << "Extracted:\n";
      cerr << nascent;
    }

    return;
  }


  // Front-end to atom extraction.  This is used to provide a
  // post-condition where an empty vGroup is returned if not all
  // torsions could be calculated for the ith residue (and the user
  // wants us to skip it)...
  vGroup extractAtoms(vGroup& residues, const int i) {
    vGroup result = extractImpl(residues, i);

    if (skip_when_missing) {
      vGroup::iterator vi;
      for (vi = result.begin(); vi != result.end(); ++vi)
        if ((*vi).size() != 4)
          return(vGroup());
    }

    return(result);
  }


protected:


  // This is used by the derived classes to make sure that going out
  // of bounds on the vector index isn't fatal...

  AtomicGroup& getResidue(vGroup& residues, const int i) {

    if (i < 0 || i >= static_cast<int>(residues.size())) {
      static AtomicGroup null_group;
      return(null_group);   // Yes, returning a handle to internal
                            // data is generally a bad idea, but we
                            // really need to return a ref, alors...
    }

    return(residues[i]);
  }

private:

  // Polymorphic extractAtoms function...  This is what changes to
  // provide more types of torsions to extract.
  virtual vGroup extractImpl(vGroup&, const int) =0;


  bool missing_atoms_warn;
  bool skip_when_missing;
};



// Class for extracting phi-psi backbone torsions for proteins...

class PhiPsi : public Extractor {
public:

  virtual ~PhiPsi() { }

  // Names, i.e. phi and psi
  vector<string> names(void) {
    vector<string> s;

    s.push_back("phi");
    s.push_back("psi");
    return(s);
  }

  
  // Implementation of function to extract the appropriate atoms for
  // phi & psi...

  vGroup extractImpl(vGroup& residues, const int i) {
    
    // Select specific atom types
    AtomNameSelector carbon("C");
    AtomNameSelector nitrogen("N");
    AtomNameSelector calpha("CA");

    // Pull out the atoms we'll use.  Since AtomicGroup::select() only
    // returns an AtomicGroup, we use this even though we're actually
    // dealing with a single atom in the group...
    AtomicGroup C_prev = getResidue(residues, i-1).select(carbon);
    AtomicGroup N = getResidue(residues, i).select(nitrogen);
    AtomicGroup CA = getResidue(residues, i).select(calpha);
    AtomicGroup C = getResidue(residues, i).select(carbon);
    AtomicGroup N_succ = getResidue(residues, i+1).select(nitrogen);
    
    // Now build an AtomicGroup corresponding to the atoms used to
    // calculate phi
    AtomicGroup phi = C_prev;
    phi.append(N);
    phi.append(CA);
    phi.append(C);
    checkAtoms("phi", phi, getResidue(residues, i));
    
    AtomicGroup psi = N;
    psi.append(CA);
    psi.append(C);
    psi.append(N_succ);
    checkAtoms("psi", psi, getResidue(residues, i));

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
  virtual ~PseudoTorsions() { }

  vector<string> names(void) {
    vector<string> s;

    s.push_back("eta");
    s.push_back("theta");
    return(s);
  }



  vGroup extractImpl(vGroup& residues, const int i) {

    // Select specific atom types
    AtomNameSelector C4P("C4'");
    AtomNameSelector P("P");

    // Extract all the atoms we'll use to build the psuedo-torsions
    AtomicGroup c4p_prev = getResidue(residues, i-1).select(C4P);
    AtomicGroup p = getResidue(residues, i).select(P);
    AtomicGroup c4p = getResidue(residues, i).select(C4P);
    AtomicGroup p_succ = getResidue(residues, i+1).select(P);
    AtomicGroup c4p_succ = getResidue(residues, i+1).select(C4P);

    // Assemble the group representing eta
    AtomicGroup eta = c4p_prev;
    eta.append(p);
    eta.append(c4p);
    eta.append(p_succ);
    checkAtoms("eta", eta, getResidue(residues, i));

    AtomicGroup theta = p;
    theta.append(c4p);
    theta.append(p_succ);
    theta.append(c4p_succ);
    checkAtoms("theta", theta, getResidue(residues, i));

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

double missing_flag = -9999;  // Value to use when a torsion can't be
                              // computed...


void parseOptions(int argc, char *argv[]) {

  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("pseudo,p", "Use RNA pseudo-torsions")
      ("missing,m", po::value<double>(&missing_flag)->default_value(-9999), "Value to use for missing torsions")
      ("warn,w", "Warn when atoms are missing a torsion")
      ("skip,s", "Skip residues where not all torsions are available")
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
      extractor = new PhiPsi;                // Assume protein

    // Configure extractor to warn upon finding atoms missing, if the
    // user wants it...
    if (vm.count("warn"))
      extractor->missingAtomsWarn();

    if (vm.count("skip"))
      extractor->skipMissingResidues();

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


  // Data-structure here is a vector of vectors of AtomicGroups.  Each
  // AtomicGroup holds the list of atoms to use for the torsion calc.
  // Each inner vector of AtomicGroups describes the torsions to
  // calculate for each residue/nucleotide.  The outer-vector is the
  // list of all residues/nucleotides to operate over...

  vvGroup torsion_atoms;
  for (int i=0; i<static_cast<int>(residues.size()); ++i) {
    vGroup atoms = extractor->extractAtoms(residues, i);
    if (!atoms.empty())
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
        double angle = missing_flag;
        if ((*vi).size() == 4)
          angle = Math::torsion((*vi)[0], (*vi)[1], (*vi)[2], (*vi)[3]);
        cout << angle << "\t";
      }
      cout << endl;
    }
  }

}
