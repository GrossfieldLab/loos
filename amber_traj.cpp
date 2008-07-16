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




#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <vector>
#include <tr1/memory>

#include <stdio.h>
#include <string.h>
#include <assert.h>


using namespace std;

#include <amber_traj.hpp>

// Scan the trajectory file to determine frame sizes and box
void AmberTraj::init(void) {
  char buf[1024];

  ifs()->getline(buf, 1024);
  frame_offset = ifs()->tellg();
  greal x, y, z;

  for (uint i=0; i<_natoms; i++) {
    *(ifs()) >> setw(8) >> x >> setw(8) >> y >> setw(8) >> z;
    frame.push_back(GCoord(x,y,z));
  }

  unsigned long fpos = ifs()->tellg();

  // This is probably not a good way of doing this???
  ifs()->getline(buf, 1024);
  stringstream ss(buf);
  double a, b, c;
  ss >> setw(8) >> a >> setw(8) >> b >> setw(8) >> c;
  if (ss.eof()) {
    fpos = ifs()->tellg();
    periodic = true;
    box = GCoord(a, b, c);
  }

  frame_size = fpos - frame_offset;

  // Now try to count the number of frames...
  _nframes = 0;
  while (!ifs()->fail()) {
    ++_nframes;
    fpos = _nframes * frame_size + frame_offset;
    ifs()->seekg(fpos);
  }

  ifs()->clear();
  ifs()->seekg(frame_offset);

  // Punt our failure check to the end...for now...
  if (ifs()->fail())
    throw(runtime_error("Unable to divine frame information from amber trajectory"));

  // This is a little hook so if we don't re-read the first frame if
  // that's the first frame requested...
  unread = true;
}


bool AmberTraj::readFrame(void) {
  greal x, y, z;

  if (unread) {
    unread = false;
    return(true);
  }

  if (ifs()->eof())
    return(false);

  for (uint i=0; i<_natoms; i++) {
    *(ifs()) >> setw(8) >> x >> setw(8) >> y >> setw(8) >> z;
    frame[i] = GCoord(x, y, z);
  }

  if (periodic) {
    greal a, b, c;
    *(ifs()) >> setw(8) >> a >> setw(8) >> b >> setw(8) >> c;
    box = GCoord(a, b, c);
  }

  return(true);
}


bool AmberTraj::readFrame(const uint i) {

  if (i == 0 && unread) {
    return(true);
  }
  unread = false;

  unsigned long fpos = i * frame_size + frame_offset;
  if (fpos < 0 || fpos >= _nframes)
    throw(runtime_error("Error- attempting to read an invalid frame from an Amber trajectory"));


  ifs()->seekg(fpos);
  if (ifs()->fail())
    throw(runtime_error("Error- cannot seek to the requested frame in an Amber trajectory"));

  return(readFrame());
}




void AmberTraj::updateGroupCoords(AtomicGroup& g) {
  AtomicGroup::Iterator iter(g);
  pAtom pa;

  while (pa = iter()) {
    uint i = pa->id() - 1;
    if (i >= _natoms)
      throw(runtime_error("Attempting to index a nonexistent atom in AmberTraj::updateGroupCoords()"));
    pa->coords(frame[i]);
  }

  if (periodic)
    g.periodicBox(box);
}
