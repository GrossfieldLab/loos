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

#if !defined(CCPDB_HPP)
#define CCPDB_HPP


#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

#include <boost/utility.hpp>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>
#include <StreamWrapper.hpp>
#include <Trajectory.hpp>

#include <pdb.hpp>


namespace loos {

  //! Class for interpreting concatenated PDB files as a Trajectory.
  /** This class reads a concatenated PDB Trajectory file as a LOOS Trajectory.
   * Each frame of the trajectory must be separated by an "END" record.
   * Since each frame is a fully-parsed PDB object, there is quite a
   * bit of overhead involved in reading CCPDB trajectories.  In
   * addition, upon instantiation, the trajectory will be scanned for
   * "END" records to build a list of seek indices for each frame.
   *
   * It is possible to get the contained PDB object out of the CCPDB,
   * but be careful of semantics that are slightly inconsistent with the
   * rest of LOOS.  See CCPDB::currentFrame() for more details.
   */

  class CCPDB : public Trajectory {
  public:
    explicit CCPDB(const std::string& s) : Trajectory(s), _natoms(0), _nframes(0) { init(); }
    explicit CCPDB(const char *p) : Trajectory(p), _natoms(0), _nframes(0) { init(); }

    virtual void rewind(void) { ifs()->clear(); ifs()->seekg(0); }
    virtual uint nframes(void) const { return(_nframes); }
    virtual uint natoms(void) const { return(_natoms); }
    virtual std::vector<GCoord> coords(void);
    virtual void updateGroupCoords(AtomicGroup& g) { g.copyCoordinates(frame); }

    virtual void seekNextFrame(void) { }
    virtual void seekFrame(const uint);
    virtual bool parseFrame(void);

    virtual bool hasPeriodicBox(void) const { return(frame.isPeriodic()); }
    virtual GCoord periodicBox(void) const { return(frame.periodicBox()); }

    //! The timestep is currently meaningless for CCPDB's, so we return
    //! a nominal 1e-3.
    virtual float timestep(void) const { return(0.001); }

    //! Returns the current frame as a PDB object.
    /** CCPDB actually stores a PDB object inside of it that represents
     *  the currently read frame.  When you request that PDB, what you
     *  get is a shared copy with the internal one, that is, the
     *  contained Atom and PeriodicBox objects are shared.  However,
     *  when you read in a new frame, the internal PDB is swapped out
     *  with a new one.  So at this point, what you are left holding is
     *  actually a copy (equivalent to a deep copy) of the previously
     *  read frame.
     *
     *  In general, since currentFrame() is not part of the Trajectory
     *  interface, you should not use it unless you need it...
     */
    PDB currentFrame(void) const { return(frame); }

  private:
    void init(void);

  private:
    uint _natoms, _nframes;
    PDB frame;
    std::vector<long> indices;
  };


}

#endif
