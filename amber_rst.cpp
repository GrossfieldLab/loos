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



#include <amber_rst.hpp>
#include <AtomicGroup.hpp>
#include <iomanip>
#include <sstream>

namespace loos {



  bool AmberRst::parseFrame(void) {
    char buf[1024];
    greal x, y, z;

    // Skip the title...
    ifs()->getline(buf, 1024);


    // If the # of atoms in the rst file mismatch the expected count,
    // return an error...
    uint na;
    *(ifs()) >> na >> current_time;
    if (na != _natoms)
      return(false);

    // This should probably be in an initialization rather than here...???
    frame.reserve(na);

    for (uint i=0; i<_natoms && !(ifs()->eof()); ++i) {
      *(ifs()) >> std::setw(12) >> x >> std::setw(12) >> y >> std::setw(12) >> z;
      frame[i] = GCoord(x, y, z);
    }

    if (ifs()->eof())
      return(false);

    // Probe for velocities or periodic box...
    greal a, b, c;
    *(ifs()) >> std::setw(12) >> a >> std::setw(12) >> b >> std::setw(12) >> c;
    *(ifs()) >> std::setw(12) >> x >> std::setw(12) >> y >> std::setw(12) >> z;

    // Check to see if there are more numbers...  If so, implies we're
    // in velocities...
    *(ifs()) >> std::setw(12) >> x >> std::setw(12) >> y >> std::setw(12) >> z;

    if (ifs()->eof()) {
      periodic = true;
      box = GCoord(a, b, c);
      return(true);
    }

    // This means we probably have velocities, so now we have to skip
    // the appropriate # of atoms and try again for the box...

    for (uint i=3; i<_natoms && !(ifs()->eof()); ++i)
      *(ifs()) >> std::setw(12) >> x >> std::setw(12) >> y >> std::setw(12) >> z;

    // Now read the box in...
    *(ifs()) >> std::setw(12) >> a >> std::setw(12) >> b >> std::setw(12) >> c;
    *(ifs()) >> std::setw(12) >> x >> std::setw(12) >> y >> std::setw(12) >> z;

    if (ifs()->eof() || ifs()->fail())
      return(true);    // Apparently, no box info...

    periodic = true;
    box = GCoord(a, b, c);

    return(true);
  }



  void AmberRst::updateGroupCoords(AtomicGroup& g) {
    AtomicGroup::iterator gi;
    pAtom pa;


    for (gi = g.begin(); gi != g.end(); ++gi) {
      uint i = (*gi)->id() - 1;
      if (i >= _natoms)
        throw(std::runtime_error("Attempting to index a nonexistent atom in AmberTraj::updateGroupCoords()"));
      (*gi)->coords(frame[i]);
    }

    if (periodic)
      g.periodicBox(box);
  }


  void AmberRst::seekNextFrameImpl(void) {
    if (!seek_flag) {
      seek_flag = true;
      return;
    }

    throw(std::logic_error("Amber RST files cannot be seeked beyond the first frame"));
  }

  void AmberRst::seekFrameImpl(const uint i) {
    if (i != 0)
      throw(std::logic_error("Amber RST files cannot be seeked beyond the first frame"));
  }
}
