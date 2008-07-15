/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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



#include <ios>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <stdlib.h>
#include <ctype.h>

#include <boost/algorithm/string.hpp>

#include <amber.hpp>
#include <Fmt.hpp>

using namespace boost;


void Amber::verifyFormat(istream& is, const string fmt) {

  string str;
  is >> str;

  char_separator<char> sep("()");
  tokenizer tokens(str, sep);
  tokenizer::iterator toks = tokens.begin();

  ++toks;
  if (*toks != fmt)
    throw(runtime_error("Bad format spec"));

}


void Amber::parseCharges(istream& is) {
  verifyFormat(is, "5E16.8");

  uint n = atoms.size();
  greal m;

  for (uint i=0; i<n; i++) {
    if (! (is >> setw(16) >> m)) {
      if (is.fail())
	throw(runtime_error("IO error while reading amber charges"));
      else
	throw(runtime_error("Invalid conversion while reading amber charges"));
    }
    atoms[i]->charge(m);
  }
}



void Amber::parseMasses(istream& is) {
  verifyFormat(is, "5E16.8");

  uint n = atoms.size();
  greal m;

  for (uint i=0; i<n; i++) {
    if (! (is >> setw(16) >> m)) {
      if (is.fail())
	throw(runtime_error("IO error while reading amber masses"));
      else
	throw(runtime_error("Invalid conversion while reading amber masses"));
    }
    atoms[i]->mass(m);
  }
}



void Amber::parseResidueLabels(istream& is) {
  verifyFormat(is, "20a4");

  for (uint i=0; i<nres; i++) {
    string s;
    if (!(is >> setw(4) >> s)) {
      if (is.fail())
	throw(runtime_error("IO error while reading residue labels"));
      else
	throw(runtime_error("Invalid conversion while reading residue labels"));
    }
    residue_labels.push_back(s);
  }
}


void Amber::parseResiduePointers(istream& is) {
  verifyFormat(is, "10I8");

  for (uint i=0; i<nres; i++) {
    int j;
    if (!(is >> setw(8) >> j)) {
      if (is.fail())
	throw(runtime_error("IO error while reading residue pointers"));
      else
	throw(runtime_error("Invalid conversion while reading residue pointers"));
    }
    residue_pointers.push_back(j);
  }
}


void Amber::assignResidues(void) {
  if (!(residue_pointers.size() == nres && residue_labels.size() == nres))
    throw(runtime_error("Unable to assign residues."));

  int curresid = 0;
  uint i = 0;
  string curresname;

  for (i=0; i<nres-1; i++) {
    ++curresid;
    curresname = residue_labels[i];
    for (uint j = residue_pointers[i]; j<residue_pointers[i+1]; j++) {
      atoms[j-1]->resid(curresid);
      atoms[j-1]->resname(curresname);
    }
  }

  // Fix the end...
  ++curresid;
  curresname = residue_labels[i];
  for (uint j = residue_pointers[i]; j<natoms; j++) {
    atoms[j]->resid(curresid);
    atoms[j]->resname(curresname);
  }
}



void Amber::parseBonds(istream& is, const int n) {
  verifyFormat(is, "10I8");

  int i, a, b, k;

  for (i=0; i<n; i++) {
    if (!(is >> a >> b >> k)) {
      if (is.fail())
	throw(runtime_error("IO error while reading bonds"));
      else
	throw(runtime_error("Invalid conversion while reading bonds"));
    }

    a /= 3;
    b /= 3;

    atoms[a]->addBond(b+1);
  }
}


void Amber::parsePointers(istream& is) {
  verifyFormat(is, "10I8");
  uint dummy;

  is >> natoms;
  is >> dummy;
  is >> nbonh;
  is >> mbona;
  for (int i=0; i<7; i++)
    is >> dummy;
  is >> nres;

  for (int i=0; i<19; i++)
    is >> dummy;

  // Now build up the atomic-group...
  if (atoms.size() != 0)
    throw(logic_error("Internal error: trying to read in an amber parmtop into a non-empty group!"));

  for (uint i=0; i<natoms; i++) {
    pAtom pa(new Atom);
    pa->id(i+1);
    atoms.push_back(pa);
  }

}


// Simply slurp up the title (for now)
void Amber::parseTitle(istream& is) {
  verifyFormat(is, "20a4");
  char buf[1024];
  
  is.getline(buf, 1024);
}


void Amber::parseAtomNames(istream& is) {
  verifyFormat(is, "20a4");

  for (uint i=0; i<natoms; i++) {
    string s;
    is >> setw(4) >> s;
    atoms[i]->name(s);
  }

  if (is.fail())
    throw(runtime_error("IO error while reading atom names"));
}


void Amber::read(istream& is) {
  char buf[1024];

  is.getline(buf,1024);
  
  bool flag = false;
  string s;
  is >> s;


  while (!(is.eof() || is.fail())) {
    if (s != "%FLAG")
      throw(runtime_error("Parse error: " + s));

    is >> s;
    if (s == "TITLE")
      parseTitle(is);
    else if (s == "POINTERS")
      parsePointers(is);
    else if (s == "ATOM_NAME")
      parseAtomNames(is);
    else if (s == "CHARGE")
      parseCharges(is);
    else if (s == "MASS")
      parseMasses(is);
    else if (s == "RESIDUE_LABEL")
      parseResidueLabels(is);
    else if (s == "RESIDUE_POINTER")
      parseResiduePointers(is);
    else if (s == "BONDS_INC_HYDROGEN")
      parseBonds(is, nbonh);
    else if (s == "BONDS_WITHOUT_HYDROGEN")
      parseBonds(is, mbona);
    else {
      while (is >> s)
	if (s == "%FLAG")
	  break;
      flag = true;
    }

    if (!flag)
      is >> s;
    else
      flag = false;
    
  }

  assignResidues();
}
