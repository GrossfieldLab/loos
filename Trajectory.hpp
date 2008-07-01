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


class Trajectory : public boost::noncopyable {
public:
  Trajectory() { }
  Trajectory(const string& s) : ifs(s) { }
  Trajectory(const char* s) : ifs(s) { }
  Trajectory(fstream& fs) : ifs(fs) { }
  virtual ~Trajectory() { }

  virtual int natoms(void) const =0;
  virtual float timestep(void) const =0;
  virtual int nframes(void) const =0;

  virtual bool readFrame(void) =0;
  virtual bool readFrame(const unsigned int i) =0;
  virtual void rewind(void) =0;

  virtual bool hasBox(void) const =0;
  virtual GCoord periodicBox(void) const =0;

  // This may require interleaving of coords and hence be particularly
  // slow...
  virtual vector<GCoord> coords(void) =0;
  
  virtual void updateGroupCoords(AtomicGroup& g) =0;

protected:
  StreamWrapper ifs;


private:
};



#endif
