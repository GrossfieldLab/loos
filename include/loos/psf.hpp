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



#if !(defined LOOS_PSF_HPP)
#define LOOS_PSF_HPP


#include <fstream>
#include <stdexcept>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>


namespace loos {

  //! Class for reading a subset of the PSF format
  /**
   * Notes:
   *
   * PSF files don't have coordinates, so the coordinates property won't be
   * set.  As a result, you'll get an UnsetProperty exception if you try to 
   * access the coordinates.  The idea is you use a PSF with something else
   * which will supply the coordinates (usually some kind of Trajectory).
   *
   * This code will read both NAMD and CHARMM PSF files, since it doesn't use
   * the atom type information anyway.
   *
   * The code extracts the atom/residue/segment names and numbers, plus mass, 
   * partial charge, and connectivity.  Higher order connectivity (angles, 
   * dihedrals, etc) are ignored.
   *
   * Atomic numbers will be deduced from the masses.  No error is
   * generated if an atomic mass is unknown to LOOS.  In order to
   * verify that all atoms have an assigned mass, use the following,
   *\code
   * bool ok = psf.allHaveProperty(Atom::anumbit);
   *\endcode
  */
  class PSF : public AtomicGroup {
  public:
    PSF() { }
    virtual ~PSF() {}

    explicit PSF(const std::string& fname) : _max_index(0), _filename(fname) {
      std::ifstream ifs(fname.c_str());
      if (!ifs)
        throw(FileOpenError(fname));
      read(ifs);
    }

    explicit PSF(std::fstream &ifs) : _max_index(0), _filename("stream") {
      read(ifs);
    }

    static pAtomicGroup create(const std::string& fname) {
      return(pAtomicGroup(new PSF(fname)));
    }

    //! Clones an object for polymorphism (see AtomicGroup::clone() for more info)
    virtual PSF* clone(void) const;

    //! Creates a deep copy (see AtomicGroup::copy() for more info)
    PSF copy(void) const;

    void read(std::istream& is);  


  private:

    PSF(const AtomicGroup& grp) : AtomicGroup(grp) { }
    void parseAtomRecord(const std::string s);  

    uint _max_index;
    std::string _filename;
  };


}

#endif
