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

#if !defined(LOOS_TINKER_ARC_HPP)
#define LOOS_TINKER_ARC_HPP


#include <iostream>
#include <string>
#include <stdexcept>

#include <loos_defs.hpp>
#include <Trajectory.hpp>

#include <tinkerxyz.hpp>


namespace loos {

  //! Class for handling Tinker ARC files (concatenation of .xyz files)
  /** This class reads a concatenated .xyz tinker trajectory (i.e. an
   *  ARC file).  In order to determine the number of frames present,
   *  the trajectory is scanned from beginning to end upon
   *  instantiation.  A list of seek indices for each frame is also
   *  built.
   *
   *  The first frame is read in by init(), so there is no explicit
   *  readFrame() during initialization.
   *
   *  There seems to be an issue with some .ARC files where reading the
   *  end of the contained TinkerXYZ object does not put the input
   *  stream into an EOF state.  So, we can't depend on checking eof()
   *  in parseFrame() to flag when we've iterated off the end.
   *  TinkerArc therefore keeps track of what index into the Trajectory
   *  it's at and uses that to check to see if it's at the end or not.
   * 
   *  It is possible to get the contained TinkerXYZ object out of a
   *  TinkerArc, but with certain caveats.  See CCPDB::currentFrame()
   *  for more details.
   */

  class TinkerArc : public Trajectory {
  public:
    explicit TinkerArc(const std::string& s) : Trajectory(s), _natoms(0), _nframes(0),
                                               current_index(0), at_end(false) { init(); }
    explicit TinkerArc(const char *p) : Trajectory(p), _natoms(0),
                                        _nframes(0), current_index(0), at_end(false) { init(); }

    explicit TinkerArc(std::istream& is) : Trajectory(is), _natoms(0),
                                        _nframes(0), current_index(0), at_end(false) { init(); }

    virtual uint nframes(void) const { return(_nframes); }
    virtual uint natoms(void) const { return(_natoms); }
    virtual std::vector<GCoord> coords(void);

    virtual bool hasPeriodicBox(void) const { return(frame.isPeriodic()); }
    virtual GCoord periodicBox(void) const { return(frame.periodicBox()); }

    virtual float timestep(void) const { return(0.001); }

    //! Returns the contained TinkerXYZ object.
    /** See CCPDB::currentFrame() for some important notes about using
     *  this function.
     */
    TinkerXYZ currentFrame(void) const { return(frame); }

  private:
    void init(void);
    virtual void rewindImpl(void) { ifs()->clear(); ifs()->seekg(0); current_index = 0; at_end = false; }
    virtual void seekNextFrameImpl(void);
    virtual void seekFrameImpl(const uint);
    virtual bool parseFrame(void);
    virtual void updateGroupCoordsImpl(AtomicGroup& g);


  private:
    uint _natoms, _nframes;
    uint current_index;
    bool at_end;
    TinkerXYZ frame;
    std::vector<long> indices;
  };


}


#endif
