/*
  rotamer.cpp


  Computes chi-1, chi-2 angles for selected side-chains.  If the
  requested angle doesn't exist, then -9999.99 is output as a marker.
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008 Tod D. Romo
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
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <cmath>
#include <sstream>
#include <limits>

using namespace std;
using namespace loos;

const double null_value = -9999.99;

/// @cond TOOLS_INTERNAL


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tCalculate the chi-1 and chi-2 angles for selected side chains\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool will calculate the chi-1 and chi-2 angles for the selected\n"
    "side chains, writing them out as a time-series (matrix).  For residues\n"
    "where an angle does not exist, -9999.99 is written as a marker.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\trotamer model.pdb simulation.dcd 'resid >= 50 && resid <= 59' >data.asc\n"
    "This example calculates chi-1 and chi-2 for residues 50 through 59.\n"
    "\n"
    "\trotamer model.pdb simulation.dcd 'segid == \"BAR2\"' >data.asc\n"
    "This example calculates chi-1 and chi-2 for all residues in the \"BAR2\" segment.\n"
    "\n"
    "NOTES\n"
    "\tThe selection must include all atoms in each residue necessary for determining\n"
    "chi-1 and chi-2.\n"
    "\n"
    "SEE ALSO\n"
    "\tramachandran\n";

  return(msg);
}


// Interface for torsion calculation class.

class Torsion {
public:
  virtual double torsion(void) =0;
  virtual ~Torsion() { }
};


// This class records specific pAtoms and will calculate their torsion
// angle...

class TorsionedAtoms : public Torsion {
public:
  TorsionedAtoms(pAtom a, pAtom b, pAtom c, pAtom d) {
    atoms = a + b + c + d;
  }
  virtual ~TorsionedAtoms() { }

  double torsion(void) { return(Math::torsion(atoms[0], atoms[1], atoms[2], atoms[3])); }

  friend ostream& operator<<(ostream& os, const TorsionedAtoms& t) {
    os << t.atoms;
    return(os);
  }

private:
  AtomicGroup atoms;
};



// This class just returns the null_value...  It's used when a torsion
// angle doesn't exist (i.e. chi-2 for Cys...)

class NoTorsion : public Torsion {
public:
  double torsion(void) { return(null_value); }
  virtual ~NoTorsion() { }
};


/// @endcond

// Utility function to pick a specific atom by name out of a group.
// It's overkill, but doing it this way does make it easy to note when
// a selection may be malformed and it only needs to be done once at
// startup...

pAtom pickAtom(const AtomicGroup& grp, const string& name) {
  AtomNameSelector sel(name);
  AtomicGroup pick = grp.select(sel);

  pAtom atom;
  if (pick.size() > 1)
    cerr << boost::format("WARNING - found more than one %s in :\n%s\n") % name % grp;
  if (pick.size() >= 1)
    atom = pick[0];
  else {
    cerr << "ERROR - could not find " << name << " in:\n" << grp << endl;
    exit(-10);
  }
  
  return(atom);
}


// Factory function for binding the torsion calculation to a group of
// atoms..

Torsion* torsionFactory(const AtomicGroup& grp, const string& a, const string& b, const string& c, const string& d) {
  if (a == "-" || b == "-" || c == "-" || d == "-")
    return(new NoTorsion());

  return(new TorsionedAtoms(pickAtom(grp, a), pickAtom(grp, b), pickAtom(grp, c), pickAtom(grp, d)));
}


// ---------------------------------------------------------------

// Map of residues to chi-1, chi-2 atom lists.  An atom name of "-"
// indicates that this angle can't be calculated.  An entry of all "-"
// marks the end of the list...

string angle_mapping[][3] = {
  { "GLY", "-,-,-,-", "-,-,-,-" },
  { "ALA", "-,-,-,-", "-,-,-,-" },
  { "VAL", "N,CA,CB,CG1", "-,-,-,-" },
  { "LEU", "N,CA,CB,CG", "CA,CB,CG,CD1" },
  { "ILE", "N,CA,CB,CG1", "CA,CB,CG1,CD" },
  { "PRO", "N,CA,CB,CG", "CA,CB,CG,CD" },
  { "PHE", "N,CA,CB,CG", "CA,CB,CG,CD1" },
  { "TYR", "N,CA,CB,CG", "CA,CB,CG,CD1" },
  { "TRP", "N,CA,CB,CG", "CA,CB,CG,CD1" },
  { "SER", "N,CA,CB,OG", "-,-,-,-" },
  { "THR", "N,CA,CB,OG1", "-,-,-,-" },
  { "CYS", "N,CA,CB,SG", "-,-,-,-" },
  { "MET", "N,CA,CB,CG", "CA,CB,CG,SD" },
  { "MSE", "N,CA,CB,CG", "CA,CB,CG,SE" },
  { "LYS", "N,CA,CB,CG", "CA,CB,CG,CD" },
  { "HIS", "N,CA,CB,CG", "CA,CB,CG,ND1" },
  { "ARG", "N,CA,CB,CG", "CA,CB,CG,CD" },
  { "ASP", "N,CA,CB,CG", "CA,CB,CG,OD1" },
  { "ASN", "N,CA,CB,CG", "CA,CB,CG,OD1" },
  { "GLN", "N,CA,CB,CG", "CA,CB,CG,CD" },
  { "GLU", "N,CA,CB,CG", "CA,CB,CG,CD" },
  { "---", "-,-,-,-", "-,-,-,-" }
};


// @cond TOOLS_INTERNAL

// Convenience class for grouping 4 strings together...

struct DihedralAtoms {
  DihedralAtoms() : a("-"), b("-"), c("-"), d("-") { }
  DihedralAtoms(const string& sa, const string& sb, const string& sc, const string& sd) :
    a(sa), b(sb), c(sc), d(sd) { }

  DihedralAtoms(const DihedralAtoms& t) : a(t.a), b(t.b), c(t.c), d(t.d) { }

  friend ostream& operator<<(ostream& os, const DihedralAtoms& t) {
    os << boost::format("(%s,%s,%s,%s)") % t.a % t.b % t.c % t.d;
    return(os);
  }

  string a, b, c, d;
};



// Mapping of residue names to the appropriate atoms for calculating
// torsion angles...

typedef boost::unordered_map<string, DihedralAtoms> ResidueDihedralAtoms;

ResidueDihedralAtoms Chi1Atoms;
ResidueDihedralAtoms Chi2Atoms;

// @endcond

// Initialize the maps...

void makeMaps(void) {

  for (int i=0; angle_mapping[i][0] != "---"; ++i) {
    vector<string> tokens;
    string name = angle_mapping[i][0];
    tokens = boost::split(tokens, angle_mapping[i][1], boost::is_any_of(","));
    if (tokens.size() != 4)
      throw(logic_error("Invalid count of dihedral atoms"));
    Chi1Atoms[name] = DihedralAtoms(tokens[0], tokens[1], tokens[2], tokens[3]);
    tokens.clear();

    tokens = boost::split(tokens, angle_mapping[i][2], boost::is_any_of(","));
    if (tokens.size() != 4)
      throw(logic_error("Invalid count of dihedral atoms"));
    Chi2Atoms[name] = DihedralAtoms(tokens[0], tokens[1], tokens[2], tokens[3]);
  }
}


// Given a map of residue to torsion atoms, pull them out of the
// passed group and bind it to a torsion calculator...

Torsion* makeTorsion(const AtomicGroup& grp, ResidueDihedralAtoms& binding) {

  // Note: G++ < 4.1 has a bug in tr1::unordered_map where there is no
  // default constructor for the iterators, so we must use the CC to
  // instantiate them...  This also apparently causes an issue when
  // trying to make binding const (and use a const_iterator)
  ResidueDihedralAtoms::iterator i(binding.end());

  string name = grp[0]->resname();
  i = binding.find(name);
  if (i == binding.end()) {
    cerr << "ERROR - no torsion information available for " << name << endl;
    exit(-20);
  }

  DihedralAtoms atoms = i->second;
  return(torsionFactory(grp, atoms.a, atoms.b, atoms.c, atoms.d));
}


// ---------------------------------------------------------------


int main(int argc, char *argv[]) {
  
  if (argc < 4) {
    cerr << "Usage - rotamer model traj sel-1 [sel-2 ...] >output.asc\n";
    cerr << fullHelpMessage();
    exit(-1);
  }

  ostringstream sshdr;
  sshdr << invocationHeader(argc, argv) << endl;
  makeMaps();

  AtomicGroup model = loos::createSystem(argv[1]);
  pTraj traj = loos::createTrajectory(argv[2], model);

  // Build the list of atoms/torsion angles to calculate...
  vector<Torsion*> chi1;
  vector<Torsion*> chi2;
  uint idx = 2;
  for (int i=3; i<argc; i++) {
    AtomicGroup subset = loos::selectAtoms(model, argv[i]);
    vector<AtomicGroup> residues = subset.splitByResidue();
    for (uint j=0; j<residues.size(); ++j, idx += 2) {
      sshdr << boost::format("# %d = %d %s %s %d %s%s")
        % idx
        % residues[j][0]->id()
        % residues[j][0]->name()
        % residues[j][0]->resname()
        % residues[j][0]->resid()
        % residues[j][0]->segid()
        % (i == argc-1 && j == residues.size()-1 ? "" : "\n");

      Torsion *x1 = makeTorsion(residues[j], Chi1Atoms);
      Torsion *x2 = makeTorsion(residues[j], Chi2Atoms);
    
      chi1.push_back(x1);
      chi2.push_back(x2);
    }
  }


  uint rows = traj->nframes();
  uint n = chi1.size();
  uint cols = 2 * n;
  Math::Matrix<double, Math::RowMajor> M(rows, cols+1);

  for (uint j=0; j<rows; ++j) {
    traj->readFrame(j);
    traj->updateGroupCoords(model);

    M(j, 0) = j;
    for (uint i = 0; i<n; i++) {
      M(j, 2*i+1) = chi1[i]->torsion();
      M(j, 2*i+2) = chi2[i]->torsion();
    }    
  }

  string hdr = sshdr.str();
  writeAsciiMatrix(cout, M, hdr);
}
