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


#include <ccpdb.hpp>

void CCPDB::init(void) {
     char buf[512];

     // Read the first frame to get the # of atoms...
     frame.read(*(ifs()));
     _natoms = frame.size();
     cached_first = true;
     indices.push_back(ifs()->tellg());

     // Now count the # of END statements...
     while (ifs()->getline(buf, sizeof(buf)))
       if (strncmp("END", buf, 3) == 0)
	 indices.push_back(ifs()->tellg());

     _nframes = indices.size();
     ifs()->clear();
     ifs()->seekg(indices[0]);
}


void CCPDB::seekFrame(const uint i) {
  if (i >= _nframes)
    throw(runtime_error("Error- Attempting to access more frames than are in the trajectory."));

  ifs()->seekg(indices[i]);
  if (ifs()->fail())
    throw(runtime_error("Error- cannot seek to the requested frame in trajectory."));
}


bool CCPDB::parseFrame(void) {
  if (ifs()->eof())
    return(false);

  // We cheat here...  Maybe it would be better to have a PDB::clear()
  // or AtomicGroup::clear() member function???
  PDB newframe;
  newframe.read(*(ifs()));
  frame = newframe;
  if (frame.size() == 0)
    return(false);

  return(true);
}


vector<GCoord> CCPDB::coords(void) {
  vector<GCoord> result(_natoms);

  for (uint i=0; i<_natoms; i++)
    result[i] = frame[i]->coords();

  return(result);
}
