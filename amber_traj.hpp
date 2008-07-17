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

#if !defined(AMBER_TRAJ_HPP)
#define AMBER_TRAJ_HPP


#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

#include <boost/utility.hpp>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>
#include <StreamWrapper.hpp>
#include <Trajectory.hpp>


//! Class for reading amber coordinate trajectories
/*!
 * This class will read in the first frame of the trajectory upon
 * instantiation.  It will also scan the file to determine how many
 * frames there are.
 *
 * Since the Amber trajectory format does not store the # of atoms
 * present, this must be passed to the AmberTraj constructor.
 *
 * Note that the Amber timestep is (presumably) defined in the parmtop
 * file, not in the trajectory file.  So we return a null-value here...
 */

class AmberTraj : public Trajectory {
public:
  explicit AmberTraj(const string& s, const int na) : Trajectory(s), _natoms(na), frame_offset(0), frame_size(0), periodic(false), unread(false) { init(); }
  explicit AmberTraj(const char* p, const int na) : Trajectory(p), _natoms(na), frame_offset(0), frame_size(0), periodic(false), unread(false) { init(); }

  virtual void rewind(void) { ifs()->seekg(frame_offset); }
  virtual uint nframes(void) const { return(_nframes); }
  virtual uint natoms(void) const { return(_natoms); }
  virtual vector<GCoord> coords(void) { return(frame); }
  virtual void updateGroupCoords(AtomicGroup&);

  virtual bool readFrame(const uint);
  //! Trajectory frame iterator
  /*!
   * After an EOF has been reached and readFrame() returns a false,
   * the cached frame is likely invalid.
   */
  virtual bool readFrame(void);

  virtual bool hasPeriodicBox(void) const { return(periodic); }
  virtual GCoord periodicBox(void) const { return(box); }

  /*!
   * As stated above, Amber does not store the timestep in the
   * trajectory, but in the parmtop instead.  So we return a
   * null-value here...
   */
  virtual float timestep(void) const { return(0.0); }  // Dummy routine...

private:
  void init(void);

private:
  uint _natoms, _nframes;
  unsigned long frame_offset, frame_size;
  bool periodic;
  bool unread;
  GCoord box;
  vector<GCoord> frame;

};

#endif
