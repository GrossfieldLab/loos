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

#if !defined(PDBTRAJ_HPP)
#define PDBTRAJ_HPP


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>

#include <boost/utility.hpp>
#include <boost/format.hpp>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>
#include <StreamWrapper.hpp>
#include <Trajectory.hpp>

#include <pdb.hpp>


class PDBTraj : public Trajectory {
public:
  explicit PDBTraj(const string& s, uint st, uint en, uint str=1) : Trajectory(), pattern(s), start(st), end(en), stride(str), _natoms(0), _nframes(0), current_index(0), at_end(false) { init(); }
  explicit PDBTraj(const char *p, uint st, uint en, uint str=1) : Trajectory(), pattern(string(p)), start(st), end(en), stride(str), _natoms(0), _nframes(0), current_index(0), at_end(false) { init(); }


  virtual void rewind(void) { seekFrame(0); }
  virtual uint nframes(void) const { return(_nframes); }
  virtual uint natoms(void) const { return(_natoms); }
  virtual vector<GCoord> coords(void);
  virtual void updateGroupCoords(AtomicGroup& g) { g.copyCoordinates(frame); }

  virtual void seekNextFrame(void);
  virtual void seekFrame(const uint);
  virtual bool parseFrame(void);

  virtual bool hasPeriodicBox(void) const { return(frame.isPeriodic()); }
  virtual GCoord periodicBox(void) const { return(frame.periodicBox()); }

  virtual float timestep(void) const { return(0.001); }

  string currentName(void) const { return(current_name); }

  PDB currentFrame(void) const { return(frame); }

private:
  void init(void);

private:
  string pattern;
  uint start, end, stride;
  uint _natoms, _nframes;
  uint current_index;
  bool at_end;
  string current_name;
  PDB frame;
};




#endif
