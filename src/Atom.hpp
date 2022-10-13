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


#if !defined(LOOS_ATOM_HPP)
#define LOOS_ATOM_HPP

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <functional>

#include <loos_defs.hpp>
#include <exceptions.hpp>
#include <Coord.hpp>

namespace loos {

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

    //! DEPRECATED exception class...use loos::UnsetProperty instead
    class UnsetProperty : public LOOSError {
    public:
      UnsetProperty() : LOOSError("Attempting to access an unset atom property") {}
      UnsetProperty(const std::string& p) : LOOSError(p) {}
      UnsetProperty(const Atom& a, const std::string& p) : LOOSError(a, p) {}
    };



    // Be careful that we don't overflow an unsigned long!
    //! Bits in the bitmask that flag what properties have actually been set.
    enum bits {
      nullbit = 0,
      coordsbit = 1,
      bondsbit = coordsbit << 1,
      massbit = bondsbit << 1,
      chargebit = massbit << 1,
      anumbit = chargebit << 1,
      flagbit = anumbit << 1,
      usr1bit = flagbit << 1,
      usr2bit = usr1bit << 1,
      usr3bit = usr2bit << 1,
      indexbit = usr3bit << 1,
      velbit = indexbit << 1
    };

    Atom() { init(); }

    //! Constructs an atom with the atomid i, atomname s, and coordinates c.
    /**
     * Constructs a new atom.
     * \param i atom-id
     * \param s atom-name
     * \param c Coordinates
     */

    Atom(const int i, const std::string s, const GCoord& c) {
      init();
      _index = 0;
      _id = i;
      _name = s;
      _coords = c;
    }


    ~Atom() { }

    // Accessors...
    int id(void) const;
    void id(const int);

    uint index(void) const;
    void index(const uint i);
  
    int resid(void) const;
    void resid(const int);

    int atomic_number(void) const;
    void atomic_number(const int);

    std::string name(void) const;
    void name(const std::string);

    std::string altLoc(void) const;
    void altLoc(const std::string);

    std::string chainId(void) const;
    void chainId(const std::string);

    std::string resname(void) const;
    void resname(const std::string);

    std::string segid(void) const;
    void segid(const std::string);

    std::string iCode(void) const;
    void iCode(const std::string);

    std::string PDBelement(void) const;
    void PDBelement(const std::string);


#if !defined(SWIG)
    //! Returns a const ref to internally stored coordinates.
    //! This returns a const ref mainly for efficiency, rather than
    //! copying the coords...
    const GCoord& coords(void) const;

    //! Returns a writable ref to the internally stored coords.
    /** This can cause problems since we track whether the coords are
     * set or not via the bitmask.  We assume that if you're accessing
     * this as non-const, your intention is to set the coords, so the
     * bit flagging coords is set automatically.
     */
    GCoord& coords(void);

    //! Returns ref to the internally stored velocities (a GCoord)
    const GCoord& velocities(void) const;
    GCoord& velocities();

#else // !defined(SWIG)

    // For python, make sure to return a copy (not a ref), otherwise we
    // get memory errors...
    GCoord coords(void) { return(_coords); }
    GCoord velocities() { return(_velocities); }

#endif // !defined(SWIG)

    //! Sets the coords to \a c
    void coords(const GCoord&);

    //! Sets the velocities
    void velocities(const GCoord&);

    double bfactor(void) const;
    void bfactor(const double);

    double occupancy(void) const;
    void occupancy(const double);

    // Note: swig requires explicit namespace on exception objects
    double charge(void) const;

    //! Sets the charge of the atom as a double.  This is NOT the PDB spec...

    void charge(const double);

    double mass(void) const;
    void mass(const double);

    int atomType() const;
    void atomType(const int);

    //! Recordname imported from the PDB for this Atom
    //! This is mainly for atoms that come from a PDB, i.e. whether or
    //! not they were an ATOM or a HETATM
    std::string recordName(void) const;
    void recordName(const std::string);

    //! Clear all stored bonds
    void clearBonds(void);
    //! Add a bond given a pAtom (extracting the atomid of the bond)
    void addBond(const pAtom&);
    //! Add a bond to an atom-id
    void addBond(const int);

    //! Deletes the specified bond.
    void deleteBond(const int);

    //! Deletes a bond by extracting the atom-id from the passed pAtom
    void deleteBond(const pAtom&);

    //! Returns a copy of the bond list.
    std::vector<int> getBonds(void) const;

    //! Sets the bonds list
    void setBonds(const std::vector<int>& list);

    bool hasBonds(void) const;

    //! Checks to see if this atom is bound to another atom
    bool isBoundTo(const int);

    bool isBoundTo(const pAtom&);

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
    bool checkProperty(const bits bitmask);

    
    //! Sets user-defined bits
    /** The user-available property bits are flagbit, usr1bit,
     * usr2bit, and usr3bit.  Attempting to set any other bits will
     * cause a logic_error exception to be thrown.
     *
     * This facility is useful if you want to tag atoms for later
     * processing without putting them in a separate group, for
     * example...
     */
    void setProperty(const bits bitmask);
    

    //! Clears user-defined bits...
    void clearProperty(const bits bitmask);

    //! Outputs an atom in pseudo-XML
    friend std::ostream& operator<<(std::ostream&, const Atom&);

    std::string asString() const;

  private:
    void init(void);
    //! Internal function for setting a bitflag
    void setPropertyBit(const bits);
    //! Internal function for clearing a bitflag
    void clearPropertyBit(const bits);

    void checkUserBits(const bits bitmask);

  private:
    int _id;
    uint _index;
    std::string _record, _name, _altloc, _resname, _chainid;
    int _resid;
    int _atomic_number;
    std::string _icode;
    double _b, _q, _charge, _mass;
    std::string _segid, _pdbelement;
    int _atom_type;
    GCoord _coords;
    GCoord _velocities;
    unsigned long mask;

    std::vector<int> bonds;
  };


#if !defined(SWIG)
  //! Compares two atoms based solely on name, id, resid, resname, and segid
  struct AtomEquals : public std::binary_function<pAtom, pAtom, bool> {
    bool operator()(const pAtom& a, const pAtom& b) const;
  };


  //! Compares two atoms based on name, id, resid, resname, segid, and coords
  /**
   * The default distance threshold is 1e-3 Angstroms
   */
  struct AtomCoordsEquals : public std::binary_function<pAtom, pAtom, bool> {
    AtomCoordsEquals(const double d) : threshold(d*d) { }
    AtomCoordsEquals() : threshold(1e-6) { }

    bool operator()(const pAtom& a, const pAtom& b) const;
    double threshold;
  };

#endif // !defined(SWIG)
}

#endif
