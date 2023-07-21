/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2023, Tod D. Romo, Alan Grossfield
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



#if !(defined LOOS_MDTRAJ_HPP)
#define LOOS_MDTRAJ_HPP

#include <ios>
#include <fstream>
#include <sstream>
#include <iostream>

#include <stdexcept>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>
#include <H5Cpp.h>
#include <boost/json.hpp>


namespace loos {

  //! Class for reading a MDTraj HDF5 file
  /**
   *
   *  This class treats the HDF5 file as a system file,
   *  and reads the topology plus the coordinates of the first frame
   *
  */
  class MDTraj : public AtomicGroup {
  public:
    MDTraj() : _max_index(0) { }
    virtual ~MDTraj() {}

    explicit MDTraj(const std::string fname) : _max_index(0), _filename(fname) {
      read();
    }

    explicit MDTraj(std::istream &ifs) : _max_index(0), _filename("stream") {
      throw(LOOSError("Creating an MDTraj from a stream isn't implemented"));
    }

    static pAtomicGroup create(const std::string& fname) {
      return(pAtomicGroup(new MDTraj(fname)));
    }

    //! Clones an object for polymorphism (see AtomicGroup::clone() for more info)
    virtual MDTraj* clone(void) const;

    //! Creates a deep copy (see AtomicGroup::copy() for more info)
    MDTraj copy(void) const;

    void read();  


  private:

    MDTraj(const AtomicGroup& grp) : AtomicGroup(grp) { }

    uint _max_index;
    std::string _filename;

    std::string getTopology(H5::H5File &file);
    void topologyToAtoms(const boost::json::value& topology);
    void topologyToBonds(const boost::json::value& topology);

  };


}

#endif
