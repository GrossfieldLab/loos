/*
 *   psf.hpp
 *   (c) 2008 Alan Grossfield and Tod D. Romo
 *   Department of Biochemistry and Biophysics
 *   University of Rochester Medical School
 *
 *   Simple CHARMM/NAMD PSF file reader
 *
 */

#if !(defined PSF_HPP)
#define PSF_HPP

#include <ios>
#include <fstream>
#include <sstream>
#include <iostream>

#include "loos.hpp"
#include "Atom.hpp"
#include "AtomicGroup.hpp"

using namespace std;

//! Class for reading a subset of the PSF format
/**
 * Notes:
 *
 *  - Coords are initialized to (99999.99, 99999.99, 99999.99) by
      default since the PSF does not contain coordinate information.
      This will hopefully make it obvious when a PSF is used without a
      matching PDB or DCD...
 */
class PSF : public AtomicGroup {
public:
    PSF() { }
    virtual ~PSF() {}

    PSF(const string fname) {
        ifstream ifs(fname.c_str());
        if (!ifs) {
            throw(runtime_error("Cannot open PSF file " + string(fname)));
            }
        read(ifs);
    }

    PSF(ifstream &ifs) {
        read(ifs);
    }

  //! Clones an object for polymorphism (see AtomicGroup::clone() for more info)
  virtual PSF* clone(void) const {
    return(new PSF(*this));
  }

  //! Creates a deep copy (see AtomicGroup::copy() for more info)
  PSF copy(void) const {
    AtomicGroup grp = this->AtomicGroup::copy();
    PSF p(grp);

    // Add PSF specific member data copies here...
    return(p);
  }

   void read(istream& is);  


private:

  PSF(const AtomicGroup& grp) : AtomicGroup(grp) { }


  void parseAtomRecord(const string s);  

};




#endif
