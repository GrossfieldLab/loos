/*
  pdb.h
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School


  Simple PDB reader/writer
*/


#if !defined(PDB_HPP)
#define PDB_HPP




#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <tr1/memory>

#include <boost/tuple/tuple.hpp>

using namespace std;
using namespace tr1;


#include <loos.hpp>
#include <Atom.hpp>
#include <AtomicGroup.hpp>
#include <pdb_remarks.hpp>

/*
  Constructors:

  o Empty
  o Reading from a filename
  o Reading from a stream
  o Created (light copy) from an AtomicGroup
*/

class PDB : public AtomicGroup {
public:
  PDB() : _show_charge(false), _auto_ter(true) { }
  virtual ~PDB() {}

  PDB(const string fname) : _show_charge(false), _auto_ter(true) {
    ifstream ifs(fname.c_str());
    if (!ifs)
      throw(runtime_error("Cannot open PDB file " + fname));
    read(ifs);
  }

  PDB(const char* fname) : _show_charge(false), _auto_ter(true) {
    ifstream ifs(fname);
    if (!ifs)
      throw(runtime_error("Cannot open PDB file " + string(fname)));
    read(ifs);
  }

  PDB(ifstream& ifs) : _show_charge(false), _auto_ter(true) { read(ifs); }

  PDB(const AtomicGroup& grp) : AtomicGroup(grp), _show_charge(false), _auto_ter(true) { }

  // Creates a deep copy...
  virtual PDB* clone(void) const {
    return(new PDB(*(this->AtomicGroup::clone())));
  }

  // Since charge is handled oddly for PDBs, we don't necessarily
  // write it to output by default...
  bool showCharge(void) const { return(_show_charge); }
  void showCharge(bool b = true) { _show_charge = b; }

  // Controls whether or not a "TER" record is automatically added
  // when printing out the group/PDB
  bool autoTerminate(void) const { return(_auto_ter); }
  void autoTerminate(bool b = true) { _auto_ter = b; }

  // Accessor for the remarks object...
  Remarks& remarks(void) { return(_remarks); }
  void remarks(const Remarks& r) { _remarks = r; }

  friend ostream& operator<<(ostream& os, const PDB& p);

  void read(istream& is);

private:

  // Internal routines for parsing...
  greal parseFloat(const string s, const unsigned int offset, const unsigned int len);
  gint parseInt(const string s, const unsigned int offset, const unsigned int len);
  string parseString(const string s, const unsigned int offset, const unsigned int len);

  int stringInt(const string s);
  
  // These will modify the PDB upon a successful parse...
  void parseRemark(const string s);
  void parseAtomRecord(const string s);
  void parseConectRecord(const string s);

  // Convert an Atom to a string representation in PDB format...
  string atomAsString(const pAtom p) const;

private:
  bool _show_charge;
  bool _auto_ter;
  Remarks _remarks;

};



#endif

