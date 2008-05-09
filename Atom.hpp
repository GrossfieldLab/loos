/*
  Atom.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Basic atom class
*/




#if !defined(ATOM_HPP)
#define ATOM_HPP


#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <tr1/memory>


using namespace std;
using namespace tr1;

#include <loos.hpp>

class Atom;
typedef shared_ptr<Atom> pAtom;


class Atom {
public:
  Atom() { init(); }
  Atom(const int i, const string s, const GCoord& c) {
    init();
    _id = i;
    _name = s;
    _coords = c;
  }

  // Accessors...
  int id(void) const { return(_id); }
  void id(const int i) { _id = i; }
  
  int resid(void) const { return(_resid); }
  void resid(const int i) { _resid = i; }

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

  const GCoord& coords(void) const { return(_coords); }
  void coords(const GCoord& c) { _coords = c; }

  double bfactor(void) const { return(_b); }
  void bfactor(const double d) { _b = d; }

  double occupancy(void) const { return(_q); }
  void occupancy(const double d) { _q = d ; }

  double charge(void) const { return(_charge); }
  void charge(const double d) { _charge = d ; }

  double mass(void) const { return(_mass); }
  void mass(const double d) { _mass = d ; }

  string recordName(void) const { return(_record); }
  void recordName(const string s) { _record = s; }

  void clearBonds(void) { bonds.clear(); }
  void addBond(const pAtom& p) { bonds.push_back(p->id()); }
  void addBond(const int i) { bonds.push_back(i); }

  void deleteBond(const int b) {
    vector<int>::iterator i = find(bonds.begin(), bonds.end(), b);
    if (i == bonds.end())
      throw(runtime_error("Attempting to delete a non-existent bond"));
    bonds.erase(i);
  }

  void deleteBond(const pAtom& p) { deleteBond(p->id()); }


  vector<int> getBonds(void) const { return(bonds); }

  bool hasBonds(void) const { return(bonds.size() != 0); }


  friend ostream& operator<<(ostream& os, const Atom& a) {
    os << "<ATOM ID='" << a._id << "' NAME='" << a._name << "' ";
    os << "RESID='" << a._resid << "' RESNAME='" << a._resname << "' ";
    os << "COORDS='" << a._coords << "' ";
    os << "ALTLOC='" << a._altloc << "' CHAINID='" << a._chainid << "' ICODE='" << a._icode << "' SEGID='" << a._segid << "' ";
    os << "B='" << a._b << "' Q='" << a._q << "' CHARGE='" << a._charge << "' MASS='" << a._mass << "'";
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
  }



private:
  int _id;
  string _record, _name, _altloc, _resname, _chainid;
  int _resid;
  string _icode;
  double _b, _q, _charge, _mass;
  string _segid, _pdbelement;
  GCoord _coords;

  vector<int> bonds;
};


#endif
