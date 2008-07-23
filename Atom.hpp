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


#if !defined(ATOM_HPP)
#define ATOM_HPP

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>


#include <loos_defs.hpp>

using namespace std;

class Atom;

//! Shared pointer to an Atom
typedef boost::shared_ptr<Atom> pAtom;

//! Basic Atom class for handling atom properties.
/**
 * This class handles atoms and atom properties.  It stores a GCoord
 * coordinate internally.  Bonds are included, but are represented as a
 * vector of atom-id's, which are assumed to be unique per atom...
 *
 * Most properties are derived from the PDB file specification.
 * Exceptions are noted below.  Accessors for each property are
 * provided and should be self-explanatory...
*/

  
class Atom {
public:

  // Be careful that we don't overflow an unsigned long!
  //! Bits in the bitmask that flag what properties have actually been set.
  enum bits {
    nullbit = 0,
    coordsbit = 1,
    bondsbit = coordsbit << 1,
    massbit = bondsbit << 1,
    chargebit = massbit << 1,
    anumbit = chargebit << 1
  };

  // Exception classes
  struct UnsetProperty : public exception {
    virtual const char *what() const throw() {
      return("Attempting to access an unset atom property.");
    }
  };

	 

  Atom() { init(); }

  //! Constructs an atom with the atomid i, atomname s, and coordinates c.
  /**
   * Constructs a new atom.
   * \param i atom-id
   * \param s atom-name
   * \param c Coordinates
  */

  Atom(const int i, const string s, const GCoord& c) {
    init();
    _id = i;
    _name = s;
    _coords = c;
  }


  ~Atom() { }

  // Accessors...
  int id(void) const { return(_id); }
  void id(const int i) { _id = i; }
  
  int resid(void) const { return(_resid); }
  void resid(const int i) { _resid = i; }

  int atomic_number(void) const { return(_atomic_number); }
  void atomic_number(const int i) { 
      _atomic_number = i;  
      setPropertyBit(anumbit);
  }

  string name(void) const { return(_name); }
  void name(const string s) { _name = s; }

  string altLoc(void) const { return(_altloc); }
  void altLoc(const string s) { _altloc = s; }

  string chainId(void) const { return(_chainid); }
  void chainId(const string s) { _chainid = s; }

  string resname(void) const { return(_resname); }
  void resname(const string s) { _resname = s; }

  string segid(void) const { return(_segid); }
  void segid(const string s) { _segid = s; }

  string iCode(void) const { return(_icode); }
  void iCode(const string s) { _icode = s; }

  string PDBelement(void) const { return(_pdbelement); }
  void PDBelement(const string s) { _pdbelement = s; }

  //! Returns a const ref to internally stored coordinates.
  //! This returns a const ref mainly for efficiency, rather than
  //! copying the coords...
  const GCoord& coords(void) const {
    if (!(mask & coordsbit))
      throw(UnsetProperty());
    return(_coords);
  }

  //! Returns a writable ref to the internally stored coords.
  /** This can cause problems since we track whether the coords are
   * set or not via the bitmask.  We assume that if you're accessing
   * this as non-const, your intention is to set the coords, so the
   * bit flagging coords is set automatically.
   */
  GCoord& coords(void) {setPropertyBit(coordsbit);  return(_coords); }

  //! Sets the coords to \a c
  void coords(const GCoord& c) { _coords = c; setPropertyBit(coordsbit); }

  double bfactor(void) const { return(_b); }
  void bfactor(const double d) { _b = d; }

  double occupancy(void) const { return(_q); }
  void occupancy(const double d) { _q = d ; }

  double charge(void) const {
    if (!(mask & chargebit))
      throw(UnsetProperty());
    return(_charge);
  }

  //! Sets the charge of the atom as a double.  This is NOT the PDB spec...

  void charge(const double d) { _charge = d ; setPropertyBit(chargebit); }

  double mass(void) const { return(_mass); }
  void mass(const double d) { _mass = d ; setPropertyBit(massbit); }

  //! Recordname imported from the PDB for this Atom
  //! This is mainly for atoms that come from a PDB, i.e. whether or
  //! not they were an ATOM or a HETATM
  string recordName(void) const { return(_record); }
  void recordName(const string s) { _record = s; }

  //! Clear all stored bonds
  void clearBonds(void) { bonds.clear(); clearPropertyBit(bondsbit); }
  //! Add a bond given a pAtom (extracting the atomid of the bond)
  void addBond(const pAtom& p) { bonds.push_back(p->id()); setPropertyBit(bondsbit); }
  //! Add a bond to an atom-id
  void addBond(const int i) { bonds.push_back(i); setPropertyBit(bondsbit); }

  //! Deletes the specified bond.
  void deleteBond(const int b) {
    vector<int>::iterator i = find(bonds.begin(), bonds.end(), b);
    if (i == bonds.end())
      throw(runtime_error("Attempting to delete a non-existent bond"));
    bonds.erase(i);
    if (bonds.size() == 0)
      clearPropertyBit(bondsbit);
  }

  //! Deletes a bond by extracting the atom-id from the passed pAtom
  void deleteBond(const pAtom& p) { deleteBond(p->id()); }

  //! Returns a copy of the bond list.
  vector<int> getBonds(void) const {
    if (!(mask & bondsbit))
      throw(UnsetProperty());
    return(bonds);
  }

  bool hasBonds(void) const { return(bonds.size() != 0); }

  //! Given a bit-mask, checks to see if those bits are set.
  /** For example, to check whether or not the coords have been set,
   *  do,
\verbatim
checkProperty(Atom::coordsbit)
\endverbatim
   *  You can combine checks by or'ing the bit flags, such as the
   *  following which checks whether both mass and charge have been
   *  set,
\verbatim
checkProperty(Atom::massbit | Atom::chargebit)
\endverbatim
  */
  bool checkProperty(const bits bitmask) { return(mask & bitmask != 0); }


  //! Outputs an atom in pseudo-XML
  friend ostream& operator<<(ostream& os, const Atom& a) {
    os << "<ATOM ID='" << a._id << "' NAME='" << a._name << "' ";
    os << "RESID='" << a._resid << "' RESNAME='" << a._resname << "' ";
    os << "COORDS='" << a._coords << "' ";
    os << "ALTLOC='" << a._altloc << "' CHAINID='" << a._chainid << "' ICODE='" << a._icode << "' SEGID='" << a._segid << "' ";
    os << "B='" << a._b << "' Q='" << a._q << "' CHARGE='" << a._charge << "' MASS='" << a._mass << "'";
    os << " ATOMICNUMBER='" << a._atomic_number <<"'";
    if (a.hasBonds() > 0) {
      vector<int>::const_iterator i;
      os << ">\n";
      for (i=a.bonds.begin(); i != a.bonds.end(); i++)
	os << "  <BOND>" << *i << "</BOND>\n";
      os << "</ATOM>";
    } else
      os << "/>";
    return(os);
  }


private:
  void init() {
    _id = -1;
    _resid = -1;
    _atomic_number = -1;
    _b = _q = 0.0;
    _charge = 0.0;
    _mass = 1.0;
    _name = "    ";
    _altloc = " ";
    _resname = "   ";
    _chainid = " ";
    _segid = "    ";
    _pdbelement = "";
    _record = "ATOM";
    mask = nullbit;   // Nullbit means nothing was set...
  }

  //! Internal function for setting a bitflag
  void setPropertyBit(const bits bitmask) { mask |= bitmask; }
  //! Internal function for clearing a bitflag
  void clearPropertyBit(const bits bitmask) { mask ^= bitmask; }

private:
  int _id;
  string _record, _name, _altloc, _resname, _chainid;
  int _resid;
  int _atomic_number;
  string _icode;
  double _b, _q, _charge, _mass;
  string _segid, _pdbelement;
  GCoord _coords;
  unsigned long mask;

  vector<int> bonds;
};


#endif
