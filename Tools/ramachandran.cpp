/*
  ramachandran.cpp



  Computes backbone torsion angles for a given set of residues...
  
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
#include <cmath>
#include <sstream>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


// Convenience typedefs...

typedef vector<AtomicGroup> vGroup;
typedef vector<vGroup> vvGroup;


// @cond TOOLS_INTERNAL

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

  Extractor() : missing_atoms_warn(false), skip_when_missing(false), show_atoms(false), verbosity_(0) { }
  virtual ~Extractor() { }

  void verbosity(uint v) { verbosity_ = v; }
  uint verbosity() const { return(verbosity_); }

  // Configuration...
  // Warn if atoms are missing and a torsion cannot be calculated...
  void missingAtomsWarn() { missing_atoms_warn = true; }

  // If atoms are missing, return an empty vGroup so this residue will
  // be skipped in future calculations...
  void skipMissingResidues() { skip_when_missing = true; }

  // Show which atoms are being used for what torsion
  void showAtoms() { show_atoms = true; }

  // Returns the names of the torsions
  virtual vector<string> names(void) =0;

  // Utility function to make sure each AtomicGroup has only 4 atoms
  // for a torsion
  void checkAtoms(const string& s, const AtomicGroup& nascent, const AtomicGroup& residue) {
    if (show_atoms && nascent.size() == 4) {
      cerr << "Extracted residues for torsion " << s << endl;
      cerr << nascent << endl;
    }

    if (!missing_atoms_warn)
      return;

    if (nascent.size() != 4) {
      cerr << boost::format("Warning- unable to determine %s from resid %d, segid '%s'\n")
        % s
        % residue[0]->resid()
        % residue[0]->segid();

      if (verbosity_ > 0) {
        cerr << "Residue dump:\n";
        cerr << residue << endl;
        cerr << "Extracted:\n";
        cerr << nascent << endl;
      }
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
        if ((*vi).size() != 4) {
          if (show_atoms)
            cerr << "***SKIPPING PREVIOUS RESIDUE***\n";
          return(vGroup());
        }
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
  bool show_atoms;
  uint verbosity_;
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


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tComputes backbone torsion angles for a given set of residues\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tGiven a set of residues, ramachandran will compute the backbone\n"
    "phi-psi angles.  The selection should include all atoms necessary to compute\n"
    "a torsion for the region of interest, i.e. it's recommended that a range of\n"
    "residues be selected by resid's and or segid's.  Not all torsions for a \n"
    "selection can be computed.  These residues are skipped in the output.  They \n"
    "can be included by using the --skipmissing=1 flag.  In this case, the missing\n"
    "torsions are replaced with a special value (default of -9999).\n"
    "\tramachandran also includes the pseudo-torsion algorithm for RNA as \n"
    "described in Wadley, Keating, Duarte, and Pyle (2007) JMB 372:942-57.\n"
    "This mode is enabled via the --pseudo=1 option.\n"
    "\tramachandran can make print out a rough secondary structure assignment\n"
    "based on phi/psi angles.  Use the --assign=1 option to turn this on.  \n"
    "Rectangular regions in the plot that roughly correspond to the clasically \n"
    "allowed regions are used to make the assignment, following discussion in:\n"
    "Hollingsworth, S. A.; Karplus, P. A. A Fresh Look at the Ramachandran Plot \n"
    "\tand the Occurrence of Standard Structures in Proteins. \n"
    "\tBioMolecular Concepts 2010, 1 (3–4), 271–283.\n"
    "\thttps://doi.org/10.1515/bmc.2010.022.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tramachandran --selection 'resid >= 1 && resid <= 100' model.psf simulation.dcd\n"
    "\n"
    "This outputs the phi-psi torsions for the first 100 residues, skipping residues\n"
    "with missing torsions.\n"
    "\n"
    "\tramachandran --selection 'resid >= 1 && resid <= 100' --assign=1 model.psf simulation.dcd\n"
    "\n"
    "This is the same as above, but also outputs the secondary structure assignment\n"
    "for each residue at each time step.\n"
    "\n"
    "\tramachandran --selection 'resid <= 200' --pseudo=1 rna.psf simulations.dcd\n"
    "\n"
    "This outputs the pseudo-torsions for the first 200 nucleic adics.\n"
    "\n"
    "\tramachandran --selection 'segid == \"PROT\"' --skipmissing=0 model.pdb simulation.dcd\n"
    "\n"
    "This outputs the phi-psi torsions for all residues in the PROT segment.  \n"
    "Residues missing torsions will have the corresponding torsion replaced \n"
    "with -9999 (default special value).\n"
    "\n"
    "NOTE: when working with proteins, as a rule we can't compute phi/psi for the \n"
    "first and last residues in your selection, because the torsions contain atoms\n" 
    "in the prior and following residues, respectively.  So, if you wanted to get \n"
    "the ramachandran map for residue 26, you'd need to use a selection like\n"
    "\n"
    "--selection '(resid >= 25) && (resid <= 27)'\n"
    "\n"
    "so that all of the atoms needed to compute the torsion for residue 26 are \n"
    "present in the selection.\n"
    ;


  return(msg);
}



// Globals...yuck!

Extractor *extractor;    // Extractor object that's used to specify
                         // which torsions to compute at run-time...

double missing_flag = -9999;  // Value to use when a torsion can't be
                              // computed...


class ToolOptions : public opts::OptionsPackage {
public:
  
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("pseudo", po::value<bool>(&pseudo_flag)->default_value(false), "Use RNA pseudo-torsions")
      ("missing", po::value<double>(&missing_flag)->default_value(-9999), "Value to use for missing torsions")
      ("warn", po::value<bool>(&warn_flag)->default_value(true), "Warn when atoms are missing a torsion")
      ("skipmissing", po::value<bool>(&skip_flag)->default_value(true), "Skip residues missing torsions")
      ("show", po::value<bool>(&show_flag)->default_value(false), "Show atoms used for each torsion")
      ("assign", po::value<bool>(&ss_flag)->default_value(false), "Assign secondary structure based on clasically allowed regions");
  }

  bool postConditions(po::variables_map& vm) {
    // Instantiate correct Extractor-derived object based on
    // user-selected mode...
    if (pseudo_flag)
      extractor = new PseudoTorsions;
    else
      extractor = new PhiPsi;                // Assume protein

    // Configure extractor to warn upon finding atoms missing, if the
    // user wants it...
    if (warn_flag)
      extractor->missingAtomsWarn();

    if (skip_flag)
      extractor->skipMissingResidues();

    if (show_flag)
      extractor->showAtoms();

    return(true);
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("pseudo=%d, missing=%f, warn=%d, skipmissing=%d, show=%d")
      % pseudo_flag % missing_flag % warn_flag % skip_flag % show_flag;
    return(oss.str());
  }

  bool pseudo_flag, warn_flag, skip_flag, show_flag;
  bool ss_flag;

};

// @endcond




// Eyeballing canonical ramachandran regions...

char classifySecondaryStructure(double phi, double psi) {

  if (phi == missing_flag || psi == missing_flag)
    return('?');

  if ( ( (psi < 0.0 && psi > -60) && phi <= -40 )
       || ( (psi > 25 && psi <= 90) && (phi >= 45 && phi <= 65) ) )
    return('H');

  if ( psi >= 90 && phi <= -45 )
    return('S');

  return('O');
}






int main(int argc, char *argv[]) {
  
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  extractor->verbosity(bopts->verbosity);

  // Read-in/setup objects to access data...
  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;
  AtomicGroup subset = selectAtoms(model, sopts->selection);

  vector<uint> indices = tropts->frameList();



  // Data-structure here is a vector of vectors of AtomicGroups.  Each
  // AtomicGroup holds the list of atoms to use for the torsion calc.
  // Each inner vector of AtomicGroups describes the torsions to
  // calculate for each residue/nucleotide.  The outer-vector is the
  // list of all residues/nucleotides to operate over...
  vvGroup torsion_atoms;

  vGroup chains = subset.splitByUniqueSegid();
  for (vGroup::iterator chain = chains.begin(); chain != chains.end(); ++chain) {

    // Split selection into individual residues to operate over...
    vGroup residues = chain->splitByResidue();

    for (int i=0; i<static_cast<int>(residues.size()); ++i) {
      vGroup atoms = extractor->extractAtoms(residues, i);
      if (!atoms.empty())
        torsion_atoms.push_back(atoms);
    }
  }

  // Verify SS option
  vector<string> torsion_names = extractor->names();
  if (torsion_names.size() != 2 ||
      !(torsion_names[0] == "phi" && torsion_names[1] == "psi")) {
    cerr << "Error- Secondary structure assignment can only be used with phi/psi torsions.\n";
    exit(-10);
  }

  
  cout << "# " << hdr << endl;
  if (topts->ss_flag)
    cout << "# Secondary Structure Codes: H = Helix, S = Sheet, O = Other, ? = Undefined\n";

  cout << "# frame\tresid" << setw(10);
  

  // Construct the header of what torsions were computed...
  copy(torsion_names.begin(), torsion_names.end(), ostream_iterator<string>(cout, "\t"));
  if (topts->ss_flag)
    cout << "SS";
  cout << endl;

  uint t = 0;
  uint resid;
  // Iterate over the requested frames from the trajectory...
  for (vector<uint>::iterator frameno = indices.begin(); frameno != indices.end(); ++frameno) {
    traj->readFrame(*frameno);

    // Update ALL atoms...since atoms are shared between AtomicGroups,
    // updating the parent group will update all of the groups we
    // built to calculate torsions...
    traj->updateGroupCoords(model);

    // Iterate over each residue...
    vvGroup::iterator vvi;
    for (vvi = torsion_atoms.begin(); vvi != torsion_atoms.end(); ++vvi) {
      cout << t << " ";

      // Iterate over each group of atoms to use for torsions within
      // the residue...
      vGroup::iterator vi;
      vector<double> torsions;

      // Grab the resid for all the torsions. Assume that the third at in first torsion is within residue.
      resid = (*(*vvi).begin())[2]->resid();
      cout << "  " << resid;
      for (vi = (*vvi).begin(); vi != (*vvi).end(); ++vi) {
        double angle = missing_flag;
        if ((*vi).size() == 4)
          angle = Math::torsion((*vi)[0], (*vi)[1], (*vi)[2], (*vi)[3]);
        cout << setw(10) << angle << "     ";
        torsions.push_back(angle);
      }
      
      if (topts->ss_flag) {
        // double-check
        if (torsions.size() != 2) {
          cerr << "Error- secondary structure requested but incorrect number of torsions found.\n";
          exit(-10);
        }
        cout << classifySecondaryStructure(torsions[0], torsions[1]);
      }
      cout << endl;
    }

    ++t;
  }

}
