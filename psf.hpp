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

  PSF(const AtomicGroup& grp) : AtomicGroup(grp) { }

    // till I figure out what's wrong with this
    virtual PSF *clone(void) const {
        return(new PSF(*(this->AtomicGroup::clone())));
    }

   void read(istream& is);  


private:
  void parseAtomRecord(const string s);  

};




#endif
