%module Atom
%include "std_string.i"
%include "std_vector.i"
%include <boost_shared_ptr.i>

%template(IntVector) std::vector<int>;
%shared_ptr(loos::Atom)

%include "Coord.i"

%{
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <functional>

#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

#include <loos_defs.hpp>
#include "Atom.hpp"

%}




namespace loos {

  class Atom;
  typedef Coord<double>   GCoord;
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
      anumbit = chargebit << 1,
      flagbit = anumbit << 1,
      usr1bit = flagbit << 1,
      usr2bit = usr1bit << 1,
      usr3bit = usr2bit << 1
    };

    Atom() { init(); }
    Atom(const int i, const std::string s, const GCoord& c) {
      init();
      _id = i;
      _name = s;
      _coords = c;
    }
    ~Atom() { }

    // Accessors...
    int id(void) const;
    void id(const int);
  
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

    const GCoord& coords(void) const;

    GCoord& coords(void);

    void coords(const GCoord&);

    double bfactor(void) const;
    void bfactor(const double);

    double occupancy(void) const;
    void occupancy(const double);

    double charge(void) const;

    void charge(const double);

    double mass(void) const;
    void mass(const double);

    std::string recordName(void) const;
    void recordName(const std::string);

    void clearBonds(void);

    void addBond(const pAtom&);

    void addBond(const int);


    void deleteBond(const int);

    void deleteBond(const pAtom&);

    std::vector<int> getBonds(void) const;

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
    //friend std::ostream& operator<<(std::ostream&, const Atom&);

  private:
    void init(void);
    //! Internal function for setting a bitflag
    void setPropertyBit(const bits);
    //! Internal function for clearing a bitflag
    void clearPropertyBit(const bits);

    void checkUserBits(const bits bitmask);

  private:
    int _id;
    std::string _record, _name, _altloc, _resname, _chainid;
    int _resid;
    int _atomic_number;
    std::string _icode;
    double _b, _q, _charge, _mass;
    std::string _segid, _pdbelement;
    GCoord _coords;
    unsigned long mask;

    std::vector<int> bonds;
  };



  // //! Compares two atoms based solely on name, id, resid, resname, and segid
  // struct AtomEquals : public std::binary_function<pAtom, pAtom, bool> {
  //   bool operator()(const pAtom& a, const pAtom& b) const;
  // };


  // //! Compares two atoms based on name, id, resid, resname, segid, and coords
  // /**
  //  * The default distance threshold is 1e-3 Angstroms
  //  */
  // struct AtomCoordsEquals : public std::binary_function<pAtom, pAtom, bool> {
  //   AtomCoordsEquals(const double d) : threshold(d*d) { }
  //   AtomCoordsEquals() : threshold(1e-6) { }

  //   bool operator()(const pAtom& a, const pAtom& b) const;
  //   double threshold;
  // };


};

