/*
  pdb.cpp
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School


  Simple PDB reader/writer
*/

#include <ios>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <stdlib.h>
#include <ctype.h>

#include <boost/algorithm/string.hpp>

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
  s << " " << setw(4) << p->name();

  s << setw(1) << p->altLoc();
  s << setw(3) << p->resname();

  s << setw(2) << p->chainId();
  s << setw(4) << p->resid();
  s << setw(2) << p->iCode();
  s << "  ";
  s << crdfmt(p->coords().x());
  s << crdfmt(p->coords().y());
  s << crdfmt(p->coords().z());
  s << bqfmt(p->occupancy());
  s << bqfmt(p->bfactor());
  s << "      ";
  s << setw(4) << p->segid();
  s << setw(2) << p->PDBelement();
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

//! Top level parser...
//! Reads a PDB from an input stream
void PDB::read(istream& is) {
  string input;

  while (getline(is, input)) {
    if (input.substr(0, 4) == "ATOM" || input.substr(0,6) == "HETATM")
      parseAtomRecord(input);
    else if (input.substr(0, 6) == "REMARK")
      parseRemark(input);
    else if (input.substr(0,6) == "CONECT")
      parseConectRecord(input);
    else if (input.substr(0,3) == "TER" || input.substr(0,3) == "END")
      ;
    else {
      int space = input.find_first_of(' ');
      cerr << "Warning- Unknown PDB record " << input.substr(0, space) << endl;
    }
  }
}


//! Output the group as a PDB...
ostream& operator<<(ostream& os, const PDB& p) {
  AtomicGroup::ConstAtomIterator i;

  os << p._remarks;
  for (i = p.atoms.begin(); i != p.atoms.end(); i++)
    os << p.atomAsString(*i) << endl;
  
  if (p._auto_ter)
    os << "TER     \n";

  return(os);
}
