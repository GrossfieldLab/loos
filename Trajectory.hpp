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

#if !defined(TRAJECTORY_HPP)
#define TRAJECTORY_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

#include <boost/utility.hpp>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>
#include <StreamWrapper.hpp>

//! Base-class for polymorphic trajectories.
/** This is the interface for trajectories in LOOS.  It is expected
 *  that at least one frame of coordinates will be buffered internally
 *  at any given time.  It is also not guaranteed that any frames will
 *  be read upon opening/instantiation.  For example, with a DCD
 *  class, the header might be read at instantiation, but no frames
 *  would necessarily be read at that point.
 *
 *  Additionally, this class is not designed to provide output, only
 *  input...
 */

class Trajectory : public boost::noncopyable {
public:
  Trajectory() { }

  //! Automatically open the file named \a s
  Trajectory(const string& s) : ifs(s) { }

  //! Automatically open the file named \a s
  Trajectory(const char* s) : ifs(s) { }

  //! Open using the given stream...
  Trajectory(fstream& fs) : ifs(fs) { }
  virtual ~Trajectory() { }

  //! # of atoms per frame
  virtual int natoms(void) const =0;
  //! Timestep per frame
  virtual float timestep(void) const =0;
  //! Number of frames in the trajectory
  virtual int nframes(void) const =0;

  //! Reading iterator
  /** This member function behaves like an iteratory.  For each call,
   * it reads the next frame into memory and returns true or, if at
   * the end of the file, returns a false.
   */
  virtual bool readFrame(void) =0;

  //! Reads a specific frame
  /** Important note: \c readFrame(i) is expected to prime the
   *iterator-like readFrame() above.
   */
  virtual bool readFrame(const unsigned int i) =0;

  //! Rewinds the readFrame() iterator
  virtual void rewind(void) =0;

  //! Tests whether or not the given frame/trajectory has periodic
  //! boundary information.
  /** The presence of periodic box information does not necessarily
   * indicate that said information has been read in yet.  For
   * example, the presence of crystal data is in the header so this
   * can be detected before any frame is read, but the crystal data
   * itself is only read when a frame is read in.
   */
  virtual bool hasPeriodicBox(void) const =0;
  //! Returns the periodic box for the current frame/trajectory
  virtual GCoord periodicBox(void) const =0;

  //! Returns the current frames coordinates as a vector of GCoords
  /** Some formats, notably DCDs, do not interleave their
   * coordinates.  This means that this could be a potentially
   * expensive operation.
   */
  virtual vector<GCoord> coords(void) =0;
  
  //! Update the coordinates in an AtomicGroup with the current frame.
  /** As with the coords() member function, some formats may have
   * non-interleaved coordinate data making the copying of coordinates
   * that much more expensive.
   *
   * In the case that the group is smaller than the trajectory frame,
   * it is assumed that the atomic-id's of the group are indices into
   * the frame's coordinates...  Since atom-id's usually begin with 1
   * and not 0, the indices are necessarily shifted by -1.
   *
   * If the trajectory has periodic boundary information, then the
   * group's periodicBox will also be updated.
   */
  virtual void updateGroupCoords(AtomicGroup& g) =0;

protected:
  StreamWrapper ifs;


private:
};



#endif
