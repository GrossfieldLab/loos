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

#if !defined(TINKER_ARC_HPP)
#define TINKER_ARC_HPP


#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

#include <boost/utility.hpp>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>
#include <StreamWrapper.hpp>
#include <Trajectory.hpp>

#include <tinkerxyz.hpp>


class TinkerArc : public Trajectory {
public:
  explicit TinkerArc(const string& s) : Trajectory(s), _natoms(0), _nframes(0) { init(); }
  explicit TinkerArc(const char *p) : Trajectory(p), _natoms(0), _nframes(0) { init(); }

  virtual void rewind(void) { ifs()->clear(); ifs()->seekg(0); }
  virtual uint nframes(void) const { return(_nframes); }
  virtual uint natoms(void) const { return(_natoms); }
  virtual vector<GCoord> coords(void);
  virtual void updateGroupCoords(AtomicGroup& g) { g.copyCoordinates(frame); }

  virtual void seekNextFrame(void) { }
  virtual void seekFrame(const uint);
  virtual bool parseFrame(void);

  virtual bool hasPeriodicBox(void) const { return(frame.isPeriodic()); }
  virtual GCoord periodicBox(void) const { return(frame.periodicBox()); }

  virtual float timestep(void) const { return(0.001); }

private:
  void init(void);

private:
  uint _natoms, _nframes;
  TinkerXYZ frame;
  vector<long> indices;
};





#endif
