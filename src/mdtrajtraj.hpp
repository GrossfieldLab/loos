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

#if !defined(LOOS_MDTRAJ_TRAJ_HPP)
#define LOOS_MDTRAJ_TRAJ_HPP


#include <string>

#include <loos_defs.hpp>
#include <Coord.hpp>
#include <Trajectory.hpp>

#include <H5Cpp.h>

namespace loos {

  //! Class for reading MDTraj HDF5 coordinate trajectories
  /*!
   * 
   */

  class MDTrajTraj : public Trajectory {
  public:
    explicit MDTrajTraj(const std::string& s, const int na) : _natoms(na), 
                                                              periodic(false),
                                                              _filename(s) { 
      init(); 
    }

    explicit MDTrajTraj(std::istream& is, const int na) : Trajectory(is), _natoms(na),
                                                          periodic(false) { 
      throw(LOOSError("Creating an MDTrajTraj from a stream isn't implemented"));
    }

    std::string description() const { return("MDTraj HDF5 trajectory"); }
    static pTraj create(const std::string& fname, const AtomicGroup& model) {
      return(pTraj(new MDTrajTraj(fname, model.size())));
    }


    virtual uint nframes(void) const { return(_nframes); }
    virtual uint natoms(void) const { return(_natoms); }
	  virtual std::vector<GCoord> coords(void) const { return(frame); }

    virtual bool hasPeriodicBox(void) const { return(periodic); }
    virtual GCoord periodicBox(void) const { return(box); }

	  //virtual double velocityConversionFactor() const { return(20.455); }

    /*!
     * MDTraj HDF5 does not store the timestep in the
     * trajectory.  So we return a
     * null-value here...
     */
    virtual float timestep(void) const { return(0.0); }  // Dummy routine...

    virtual bool parseFrame(void);


  private:
    void init(void);
    virtual void rewindImpl(void) { _current_frame = 0; }
    virtual void seekNextFrameImpl(void) { _current_frame++; seekFrameImpl(_current_frame);}
    virtual void seekFrameImpl(const uint);
    virtual void updateGroupCoordsImpl(AtomicGroup&);
    void readRawFrame(const uint i);


  private:
    uint _natoms, _nframes;
    std::string _filename;
    bool periodic;
    GCoord box;
    std::vector<GCoord> frame;

    // All of the internal stuff to read the HDF5 file
    H5::H5File file;
    H5::DataSet box_dataset;
    H5::DataSpace box_dataspace;
    H5::DataType box_datatype;
    H5::DataSpace box_memspace;
    H5::DataSet coords_dataset;
    H5::DataSpace coords_dataspace;
    H5::DataType coords_datatype;
    H5::DataSpace coord_memspace;
  };


}
#endif
