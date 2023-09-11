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



#include <amber_traj.hpp>
#include <AtomicGroup.hpp>
#include <iomanip>
#include <sstream>

namespace loos {

  // Scan the trajectory file to determine frame sizes and box
  void AmberTraj::init(void) {
    char buf[1024];

    ifs->getline(buf, 1024);
    frame_offset = ifs->tellg();
    greal x, y, z;

    for (uint i=0; i<_natoms; i++) {
      *(ifs) >> std::setw(8) >> x >> std::setw(8) >> y >> std::setw(8) >> z;
      frame.push_back(GCoord(x,y,z));
    }


    // This is probably not a good way of doing this???
    ifs->getline(buf, 1024);
    unsigned long fpos = ifs->tellg();

    ifs->getline(buf, 1024);
    if (ifs->fail())
      throw(FileOpenError(_filename, "Problem scanning Amber Trajectory"));

    std::stringstream ss(buf);
    double a= -1, b= -1, c= -1;
    ss >> std::setw(8) >> a >> std::setw(8) >> b >> std::setw(8) >> c;
    if (ss.eof()) {
      fpos = ifs->tellg();
      periodic = true;
      box = GCoord(a, b, c);
    }

    frame_size = fpos - frame_offset;

    // Now try to count the number of frames...
    _nframes = 1;
    double dummy;
    while (!ifs->fail()) {
      ++_nframes;
      fpos = _nframes * frame_size + frame_offset;
      ifs->seekg(fpos);
      *(ifs) >> dummy;
    }


    ifs->clear();
    ifs->seekg(frame_offset + frame_size);

    // Punt our failure check to the end...for now...
    if (ifs->fail())
      throw(FileOpenError(_filename, "Cannot determine frame information for Amber trajectory"));

    cached_first = true;
  }


  bool AmberTraj::parseFrame(void) {
    greal x, y, z;

    if (ifs->eof())
      return(false);

    // It seems that it's possible that, at the end of the trajectory,
    // we may have read the last datum from the last frame but will not
    // be physically at the EOF.  So, we also check inside the coord
    // read loop to see if there is an EOF.  Rather than flag this as an
    // error condition, just return a false indicating we've read past
    // the end...

    for (uint i=0; i<_natoms && !(ifs->eof()); i++) {
      *(ifs) >> std::setw(8) >> x >> std::setw(8) >> y >> std::setw(8) >> z;
      frame[i] = GCoord(x, y, z);
    }

    if (ifs->eof())
      return(false);

    if (periodic) {
      greal a, b, c;
      *(ifs) >> std::setw(8) >> a >> std::setw(8) >> b >> std::setw(8) >> c;
      box = GCoord(a, b, c);
    }

    if (ifs->fail())
      throw(FileReadError(_filename, "Problem reading from Amber trajectory"));

    return(true);
  }


  void AmberTraj::seekFrameImpl(const uint i) {

    cached_first = false;
    unsigned long fpos = i * frame_size + frame_offset;
    if (i >= _nframes)
      throw(FileError(_filename, "Attempting seek frame beyond end of trajectory"));


    ifs->seekg(fpos);
    if (ifs->fail())
      throw(FileError(_filename, "Cannot seek to frame"));
  }


  void AmberTraj::updateGroupCoordsImpl(AtomicGroup& g) {

    for (AtomicGroup::iterator i = g.begin(); i != g.end(); ++i) {
      uint idx = (*i)->index();
      if (idx >= _natoms)
        throw(LOOSError(_filename, **i, "Atom index into trajectory is out of bounds"));
      (*i)->coords(frame[idx]);
    }
    
    if (periodic)
      g.periodicBox(box);
  }
}
