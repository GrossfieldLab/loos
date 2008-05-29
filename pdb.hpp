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
#include <cryst.hpp>


//! PDB reading/writing class
class PDB : public AtomicGroup {
public:
  PDB() : _show_charge(false), _auto_ter(true), _has_cryst(false) { }
  virtual ~PDB() {}

  //! Read in PDB from a filename
  PDB(const string fname) : _show_charge(false), _auto_ter(true), _has_cryst(false) {
    ifstream ifs(fname.c_str());
    if (!ifs)
      throw(runtime_error("Cannot open PDB file " + fname));
    read(ifs);
  }

  //! Read in a PDB from a filename
  PDB(const char* fname) : _show_charge(false), _auto_ter(true), _has_cryst(false) {
    ifstream ifs(fname);
    if (!ifs)
      throw(runtime_error("Cannot open PDB file " + string(fname)));
    read(ifs);
  }

  //! Read in a PDB from an ifstream
  PDB(ifstream& ifs) : _show_charge(false), _auto_ter(true), _has_cryst(false) { read(ifs); }

  //! Create a PDB from an AtomicGroup (i.e. upcast)
  PDB(const AtomicGroup& grp) : AtomicGroup(grp), _show_charge(false), _auto_ter(true), _has_cryst(false) { }

  //! Creates a deep copy.
  virtual PDB* clone(void) const {
    return(new PDB(*(this->AtomicGroup::clone())));
  }

  bool showCharge(void) const { return(_show_charge); }
  //! Special handling for charges since the PDB form is daft
  void showCharge(bool b = true) { _show_charge = b; }

  bool autoTerminate(void) const { return(_auto_ter); }

  //! Automatically insert TER record at end of output.
  /*! Controls whether or not a "TER" record is automatically added
      when printing out the group/PDB
  */
  void autoTerminate(bool b = true) { _auto_ter = b; }

  //! Accessor for the remarks object...
  Remarks& remarks(void) { return(_remarks); }
  //! Accessor for the remarks object...
  void remarks(const Remarks& r) { _remarks = r; }

  UnitCell& unitCell(void) { return(cell); }
  void unitCell(const UnitCell& c) { cell = c; }

  //! Output as a PDB
  friend ostream& operator<<(ostream& os, const PDB& p);

  //! Read in PDB from an ifstream
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
  void parseCryst1Record(const string s);

  // Convert an Atom to a string representation in PDB format...
  string atomAsString(const pAtom p) const;

private:
  bool _show_charge;
  bool _auto_ter;
  bool _has_cryst;
  Remarks _remarks;
  UnitCell cell;

};



#endif

