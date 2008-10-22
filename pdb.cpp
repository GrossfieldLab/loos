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
#include <boost/format.hpp>

#include <pdb.hpp>
#include <Fmt.hpp>

using namespace boost;

/* The following are utilities to extract a substring and parse it as
  a specific type.  If the offset lies outside the passed string, then
  a default null value is returned.  This is in case some PDBs don't
  have a SEGID, or CHARGE, etc.
*/

greal PDB::parseFloat(const string s, const unsigned int offset, const unsigned int len) {
  greal result;

  if (offset+len > s.size())
    return(0.0);

  string t = s.substr(offset, len);
  if (!(stringstream(t) >> result))
    throw(runtime_error("Cannot parse " + t + " as a floating value"));

  return(result);
}


gint PDB::parseInt(const string s, const int unsigned offset, const unsigned int len) {
  gint result;


  if (offset+len > s.size())
    return(0);
  
  string t = s.substr(offset, len);
  if (!(stringstream(t) >> result))
    throw(runtime_error("Cannot parse " + t + " as a integer value"));

  return(result);
}



string PDB::parseString(const string s, const unsigned int offset, const unsigned int len) {

  if (offset+len > s.size())
    return(string(""));

  string t = s.substr(offset, len);
  trim(t);

  return(t);
}


// Special handling for REMARKs to ignore the line code, if
// present... 

void PDB::parseRemark(const string s) {
  string t;

  if (s[6] == ' ' && isdigit(s[7]))
    t = s.substr(11, 58);
  else
    t = s.substr(7, 62);

  _remarks.add(t);
}


// Parse an ATOM or HETATM record...

void PDB::parseAtomRecord(const string s) {
  greal r;
  gint i;
  string t;
  GCoord c;
  pAtom pa(new Atom);
  
  t = parseString(s, 0, 6);
  pa->recordName(t);

  i = parseInt(s, 6, 5);
  pa->id(i);

  t = parseString(s, 12, 4);
  pa->name(t);

  t = parseString(s, 16, 1);
  pa->altLoc(t);

  t = parseString(s, 17, 4);
  pa->resname(t);

  t = parseString(s, 21, 1);
  pa->chainId(t);

  i = parseInt(s, 22, 4);
  pa->resid(i);
  
  t = parseString(s, 26, 1);

  // Special handling of resid field since it may be frame-shifted by
  // 1 col in some cases...
  if (strictness_policy) {
    char c = t[0];
    
    if (c != ' ' && !isalpha(c))
      throw(runtime_error("Non-alpha character in iCode column of PDB"));
	
  } else {
    char c = t[0];

    if (c != ' ' && isdigit(c)) {
      i = parseInt(s, 22, 5);
      pa->resid(i);
      t = " ";
    }

  }
  pa->iCode(t);

  c[0] = parseFloat(s, 30, 8);
  c[1] = parseFloat(s, 38, 8);
  c[2] = parseFloat(s, 46, 8);
  pa->coords(c);
  
  r = parseFloat(s, 54, 6);
  pa->occupancy(r);

  r = parseFloat(s, 60, 6);
  pa->bfactor(r);

  t = parseString(s, 72, 4);
  pa->segid(t);

  t = parseString(s, 76, 2);
  pa->PDBelement(t);


  // Charge is not currently handled...
  t = parseString(s, 78, 2);

  append(pa);
}



// Convert an Atom to a string with a PDB format...

string PDB::atomAsString(const pAtom p) const {
  ostringstream s;

  // Float formatter for coords
  Fmt crdfmt(3);
  crdfmt.width(8);
  crdfmt.right();
  crdfmt.trailingZeros(true);
  crdfmt.fixed();

  // Float formatter for B's and Q's...
  Fmt bqfmt(2);
  bqfmt.width(6);
  bqfmt.right();
  bqfmt.trailingZeros(true);
  bqfmt.fixed();

  s << setw(6) << left << p->recordName();
  s << setw(5) << right << p->id();
  s << " " << setw(4) << left << p->name();

  s << setw(1) << p->altLoc();
  s << setw(4) << left << p->resname();

  s << setw(1) << right << p->chainId();
  s << setw(4) << p->resid();
  s << setw(2) << p->iCode();
  if (p->resid() < 10000)
    s << "  ";
  else
    s << " ";   // HACK!
  s << crdfmt(p->coords().x());
  s << crdfmt(p->coords().y());
  s << crdfmt(p->coords().z());
  s << bqfmt(p->occupancy());
  s << bqfmt(p->bfactor());
  s << "      ";
  s << setw(4) << left << p->segid();
  s << setw(2) << right << p->PDBelement();
  if (_show_charge)
    s << setw(2) << p->charge();
  else
    s << "  ";

  return(s.str());
}


int PDB::stringInt(const string s) {
  int i;

  if (!(stringstream(s) >> i))
    throw(runtime_error("Cannot convert " + s + " to an integer."));

  return(i);
}


// Parse CONECT records, updating the referenced atoms...
void PDB::parseConectRecord(const string s) {
  vector<string> ary;
  split(ary, s, is_any_of(" "));

  if (ary.size() < 2)
    throw(runtime_error("Cannot parse CONECT record...must have two or more entries"));

  int donor_id = stringInt(ary[0]);
  pAtom donor = findById(donor_id);
  if (donor == 0)
    throw(runtime_error("Cannot find donor atom for CONECT record"));

  vector<string>::iterator i;
  for (i = ary.begin()+1; i != ary.end(); i++) {
    int acceptor_id = stringInt(*i);
    pAtom acceptor = findById(acceptor_id);
    if (acceptor == 0)
      throw(runtime_error("Cannot find acceptor atom for CONECT record"));
    donor->addBond(acceptor);
  }
}


void PDB::parseCryst1Record(const string s) {
  greal r;
  gint i;
  string t;

  r = parseFloat(s, 6, 9);
  cell.a(r);
  r = parseFloat(s, 15, 9);
  cell.b(r);
  r = parseFloat(s, 24, 9);
  cell.c(r);

  r = parseFloat(s, 33, 7);
  cell.alpha(r);
  r = parseFloat(s, 40, 7);
  cell.beta(r);
  r = parseFloat(s, 47, 7);
  cell.gamma(r);

  // Special handling in case of mangled CRYST1 record...
  if (s.length() < 66) {
    t = s.substr(55);
    cell.spaceGroup(t);
    cell.z(-1);   // ??? 
  } else {
    t = parseString(s, 55, 11);
    cell.spaceGroup(t);
    i = parseInt(s, 66, 4);
    cell.z(i);
  }

  _has_cryst = true;
}


//! Top level parser...
//! Reads a PDB from an input stream
void PDB::read(istream& is) {
  string input;
  bool has_cryst = false;

  while (getline(is, input)) {
    if (input.substr(0, 4) == "ATOM" || input.substr(0,6) == "HETATM")
      parseAtomRecord(input);
    else if (input.substr(0, 6) == "REMARK")
      parseRemark(input);
    else if (input.substr(0,6) == "CONECT")
      parseConectRecord(input);
    else if (input.substr(0, 6) == "CRYST1") {
      parseCryst1Record(input);
      has_cryst = true;
    } else if (input.substr(0,3) == "TER")
      ;
    else if (input.substr(0,3) == "END")
      break;
    else {
      int space = input.find_first_of(' ');
      cerr << "Warning- Unknown PDB record " << input.substr(0, space) << endl;
    }
  }

  // Do some post-extraction...
  if (remarksHasBox(_remarks)) {
    GCoord c = boxFromRemarks(_remarks);
    periodicBox(c);
  } else if (has_cryst) {
    GCoord c(cell.a(), cell.b(), cell.c());
    periodicBox(c);
  }
}


ostream& FormattedUnitCell(ostream& os, const UnitCell& u) {
  os << "CRYST1";
  Fmt dists(3);
  dists.width(9).right().trailingZeros(true).fixed();
  Fmt angles(2);
  angles.width(7).right().trailingZeros(true).fixed();

  os << dists(u.a()) << dists(u.b()) << dists(u.c());
  os << angles(u.alpha()) << angles(u.beta()) << angles(u.gamma());
  os << " " << setw(10) << left << u.spaceGroup() << setw(4) << u.z();

  return(os);
}

ostream& XTALLine(ostream& os, const GCoord& box) {
    os << "REMARK  XTAL "
       << box.x() << " "
       << box.y() << " "
       << box.z();
    return(os);
}


ostream& FormatConectRecords(ostream& os, PDB& p) {
  AtomicGroup::AtomIterator ci;

  // We first have to make sure that the base AtomicGroup is sorted
  // since we will be verifying bound atoms exist by searching for
  // their ID...this would force a sort while we have iterators
  // pointing to the vector of atoms which would be bad...
  p.sort();
  
  for (ci = p.atoms.begin(); ci != p.atoms.end(); ++ci) {
    if ((*ci)->checkProperty(Atom::bondsbit)) {
      int donor = (*ci)->id();

      os << format("CONECT%5d") % donor;
      int i = 0;

      vector<int> bonds = (*ci)->getBonds();
      vector<int>::const_iterator cj;
      for (cj = bonds.begin(); cj != bonds.end(); ++cj) {
	if (++i > 4) {
	  i = 1;
	  os << format("\nCONECT%5d") % donor;
	}
	int bound_id = *cj;
	pAtom pa = p.findById(bound_id);
	if (pa == 0)
	  throw(PDB::BadConnectivity());
	os << format("%5d") % bound_id;
      }
      os << endl;
    }
  }

  return(os);
}


//! Output the group as a PDB...
ostream& operator<<(ostream& os, PDB& p) {
  AtomicGroup::AtomIterator i;

  os << p._remarks;
  if (p.isPeriodic())
    XTALLine(os, p.periodicBox()) << endl;
  if (p._has_cryst) 
    FormattedUnitCell(os, p.cell) << endl;
  for (i = p.atoms.begin(); i != p.atoms.end(); ++i)
    os << p.atomAsString(*i) << endl;

  if (p.hasBonds()) {
    int maxid = 0;
    for (i = p.atoms.begin(); i != p.atoms.end(); ++i)
      if ((*i)->id() > maxid)
	maxid = (*i)->id();

    if (maxid <= 99999)
      FormatConectRecords(os, p);
  }
  
  if (p._auto_ter)
    os << "TER     \n";

  return(os);
}
