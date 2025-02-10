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

#include <Atom.hpp>
#include <algorithm>
#include <boost/format.hpp>

namespace loos {

  int Atom::id(void) const { return(_id); }
  void Atom::id(const int i) { _id = i; }

  uint Atom::index(void) const 
  {
    return(_index);
  }
  

  void Atom::index(const uint i) 
  {
    _index = i;
    setPropertyBit(indexbit);
  }
  
  
  int Atom::resid(void) const { return(_resid); }
  void Atom::resid(const int i) { _resid = i; }

  int Atom::atomic_number(void) const { return(_atomic_number); }
  void Atom::atomic_number(const int i) { 
    _atomic_number = i;  
    setPropertyBit(anumbit);
  }

  std::string Atom::name(void) const { return(_name); }
  void Atom::name(const std::string s) { _name = s; }

  std::string Atom::altLoc(void) const { return(_altloc); }
  void Atom::altLoc(const std::string s) { _altloc = s; }

  std::string Atom::chainId(void) const { return(_chainid); }
  void Atom::chainId(const std::string s) { _chainid = s; }

  std::string Atom::resname(void) const { return(_resname); }
  void Atom::resname(const std::string s) { _resname = s; }

  std::string Atom::segid(void) const { return(_segid); }
  void Atom::segid(const std::string s) { _segid = s; }

  std::string Atom::iCode(void) const { return(_icode); }
  void Atom::iCode(const std::string s) { _icode = s; }

  std::string Atom::PDBelement(void) const { return(_pdbelement); }
  void Atom::PDBelement(const std::string s) { _pdbelement = s; }

  const GCoord& Atom::coords(void) const { return(_coords); }
  GCoord& Atom::coords(void) { setPropertyBit(coordsbit); return(_coords); }
  void Atom::coords(const GCoord& c) { _coords = c; setPropertyBit(coordsbit); }


  const GCoord& Atom::velocities() const { return(_velocities); }
  GCoord& Atom::velocities() { setPropertyBit(velbit); return(_velocities); }
  void Atom::velocities(const GCoord& v) { _velocities = v; setPropertyBit(velbit); }

  double Atom::bfactor(void) const { return(_b); }
  void Atom::bfactor(const double d) { _b = d; }

  double Atom::occupancy(void) const { return(_q); }
  void Atom::occupancy(const double d) { _q = d ; }

  double Atom::charge(void) const {
    if (!(mask & chargebit))
      throw(loos::UnsetProperty("Atom has no charge set"));
    return(_charge);
  }

    //! Sets the charge of the atom as a double.  This is NOT the PDB spec...

  void Atom::charge(const double d) { _charge = d ; setPropertyBit(chargebit); }

  double Atom::mass(void) const { return(_mass); }
  void Atom::mass(const double d) { _mass = d ; setPropertyBit(massbit); }

    //! Recordname imported from the PDB for this Atom
    //! This is mainly for atoms that come from a PDB, i.e. whether or
    //! not they were an ATOM or a HETATM
  std::string Atom::recordName(void) const { return(_record); }
  void Atom::recordName(const std::string s) { _record = s; }

    //! Clear all stored bonds
  void Atom::clearBonds(void) { bonds.clear(); clearPropertyBit(bondsbit); }
    //! Add a bond given a pAtom (extracting the atomid of the bond)
  void Atom::addBond(const pAtom& p) { bonds.push_back(p->id()); setPropertyBit(bondsbit); }
    //! Add a bond to an atom-id
  void Atom::addBond(const int i) { bonds.push_back(i); setPropertyBit(bondsbit); }

    //! Deletes the specified bond.
  void Atom::deleteBond(const int b) {
    std::vector<int>::iterator i = find(bonds.begin(), bonds.end(), b);
    if (i == bonds.end())
      throw(LOOSError(*this, "Attempting to delete a non-existent bond"));
    bonds.erase(i);
    if (bonds.size() == 0)
      clearPropertyBit(bondsbit);
  }

    //! Deletes a bond by extracting the atom-id from the passed pAtom
  void Atom::deleteBond(const pAtom& p) { deleteBond(p->id()); }

    //! Returns a copy of the bond list.
  std::vector<int> Atom::getBonds(void) const {
    if (!(mask & bondsbit))
      throw(loos::UnsetProperty("Atom has no connectivity"));
    return(bonds);
  }

  void Atom::setBonds(const std::vector<int>& list) {
    bonds = list;
    setPropertyBit(bondsbit);
  }

  bool Atom::hasBonds(void) const { return(bonds.size() != 0); }

    //! Checks to see if this atom is bound to another atom
  bool Atom::isBoundTo(const int i) {
    std::vector<int>::iterator found = find(bonds.begin(), bonds.end(), i);
    return(found != bonds.end());
  }

  bool Atom::isBoundTo(const pAtom& p) { return(isBoundTo(p->id())); }

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
  bool Atom::checkProperty(const bits bitmask) { return((mask & bitmask) != 0); }


  // DEPRECATED: will likely be removed in future versions of LOOS
  void Atom::checkUserBits(const bits bitmask) {
    if (! (bitmask & (flagbit|usr1bit|usr2bit|usr3bit)) )
      throw(LOOSError("Attempting to set a non-user property bit in an Atom"));
  }

  void Atom::setProperty(const bits bitmask) {
    setPropertyBit(bitmask);
  }

  void Atom::clearProperty(const bits bitmask) {
    clearPropertyBit(bitmask);
  }

  void Atom::atomType(const int i) {
    _atom_type = i;
  }

  int Atom::atomType() const {
    return(_atom_type);
  }


  void Atom::init() {
    _id = 1;
    _index = 0;
    _resid = 1;
    _atomic_number = -1;
    _b = 0.0;
    _q = 1.0;
    _charge = 0.0;
    _mass = 1.0;
    _name = "    ";
    _altloc = " ";
    _resname = "   ";
    _chainid = " ";
    _segid = "    ";
    _pdbelement = "";
    _record = "ATOM";
    _atom_type = -1;
    mask = nullbit;   // Nullbit means nothing was set...
  }

  void Atom::setPropertyBit(const bits bitmask) { mask |= bitmask; }
    //! Internal function for clearing a bitflag
  void Atom::clearPropertyBit(const bits bitmask) { mask &= (~bitmask); }


  std::ostream& operator<<(std::ostream& os, const loos::Atom& a) {
    os << "<ATOM INDEX='" << a._index << "' ID='" << a._id << "' NAME='" << a._name << "' ";
    os << "RESID='" << a._resid << "' RESNAME='" << a._resname << "' ";
    os << "COORDS='" << a._coords << "' ";
    os << "VELOCITIES='" << a._velocities << "' ";
    os << "ALTLOC='" << a._altloc << "' CHAINID='" << a._chainid << "' ICODE='" << a._icode << "' SEGID='" << a._segid << "' ";
    os << "B='" << a._b << "' Q='" << a._q << "' CHARGE='" << a._charge << "' MASS='" << a._mass << "'";
    os << " ATOMICNUMBER='" << a._atomic_number <<"'";
    os << " MASK='" << boost::format("%x") % a.mask << "'";
    if (a.hasBonds() > 0) {
      std::vector<int>::const_iterator i;
      os << ">\n";
      for (i=a.bonds.begin(); i != a.bonds.end(); i++)
        os << "  <BOND>" << *i << "</BOND>\n";
      os << "</ATOM>";
    } else
      os << "/>";

    return(os);
  }


  std::string Atom::asString() const {
    std::ostringstream oss;

    oss << *this;
    return oss.str();
  }
  
  bool AtomEquals::operator()(const pAtom& a, const pAtom& b) const {
    return(a->name() == b->name()
           && a->id() == b->id()
           && a->resname() == b->resname()
           && a->resid() == b->resid()
           && a->segid() == b->segid());
  }

  bool AtomCoordsEquals::operator()(const pAtom& a, const pAtom& b) const {
    bool bb = (a->name() == b->name()
               && a->id() == b->id()
               && a->resname() == b->resname()
               && a->resid() == b->resid()
               && a->segid() == b->segid());
    if (!bb)
      return(false);

    double d = a->coords().distance2(b->coords());
    return( d <= threshold );
  }

}

