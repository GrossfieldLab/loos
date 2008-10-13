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


#include <pdbtraj.hpp>

void PDBTraj::init(void) {
  seekFrame(0);
  parseFrame();
  _natoms = frame.size();
  _nframes = (end - start) / stride + 1;
}


void PDBTraj::seekFrame(const uint i) {

  uint idx = i * stride + start;
  if (idx < start || i > end)
    throw(runtime_error("Error- Attempting to access more frames than are in the trajectory."));

  stringstream s;
  s << boost::format(pattern) % idx;
  current_name = s.str();
  ifs.setStream(current_name);
  current_index = i;
  at_end = false;
}


void PDBTraj::seekNextFrame(void) {
  if (at_end)
    return;

  if (current_index >= _nframes)
    at_end = true;
  else {
    seekFrame(current_index);
    ++current_index;
  }
}


bool PDBTraj::parseFrame(void) {

  if (at_end)
    return(false);
  
  PDB newframe;
  newframe.read(*(ifs()));
  frame = newframe;
  if (frame.size() == 0) {
    at_end = true;
    return(false);
  }

  return(true);
}



vector<GCoord> PDBTraj::coords(void) {
  vector<GCoord> result(_natoms);

  for (uint i=0; i<_natoms; i++)
    result[i] = frame[i]->coords();

  return(result);
}

