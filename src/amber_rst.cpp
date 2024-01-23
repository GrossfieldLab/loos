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
    std::string buf;
    greal x, y, z;

    // Skip the title...
    std::getline(*(ifs), buf);
    if (ifs->eof())
      return(false);   //  Catch the initial EOF...

    // Check the # of atoms...
    // Note: the Amber spec says that this line should contain both
    // the number of atoms and the current time, but some restart
    // files "in the wild" only have the number of atoms...
    std::getline(*(ifs), buf);
    std::istringstream iss(buf);
    uint na;
    iss >> na;
    if (na != _natoms)
      throw(FileReadError(_filename, "Number of atoms mismatch in Amber restart file"));

    iss >> current_time;


    // This should probably be in an initialization rather than here...???
    frame.reserve(na);

    uint i;
    for (i=0; i<_natoms && !(ifs->eof()); ++i) {
      *(ifs) >> std::setw(12) >> x >> std::setw(12) >> y >> std::setw(12) >> z;
      frame[i] = GCoord(x, y, z);
    }

    if (i != _natoms)
      throw(FileReadError(_filename, "Number of atoms read is not what was expected"));

    // Probe for velocities or periodic box...
    greal a, b, c;
    *(ifs) >> std::setw(12) >> a >> std::setw(12) >> b >> std::setw(12) >> c;
    if (ifs->eof())
      return(true);

    // Finish reading the putative box...
    *(ifs) >> std::setw(12) >> x >> std::setw(12) >> y >> std::setw(12) >> z;

    // Check to see if there are more numbers...  If so, implies we're
    // in velocities...
    *(ifs) >> std::setw(12) >> x >> std::setw(12) >> y >> std::setw(12) >> z;

    if (ifs->eof()) {
      periodic = true;
      box = GCoord(a, b, c);
      return(true);
    }

    // This means we probably have velocities, so now we have to skip
    // the appropriate # of atoms and try again for the box...

    for (uint i=3; i<_natoms && !(ifs->eof()); ++i)
      *(ifs) >> std::setw(12) >> x >> std::setw(12) >> y >> std::setw(12) >> z;

    // Now read the box in...
    *(ifs) >> std::setw(12) >> a >> std::setw(12) >> b >> std::setw(12) >> c;
    *(ifs) >> std::setw(12) >> x >> std::setw(12) >> y >> std::setw(12) >> z;

    if (ifs->eof() || ifs->fail())
      return(true);

    periodic = true;
    box = GCoord(a, b, c);

    return(true);
  }



  void AmberRst::updateGroupCoordsImpl(AtomicGroup& g) {
    AtomicGroup::iterator gi;
    pAtom pa;


    for (gi = g.begin(); gi != g.end(); ++gi) {
      uint i = (*gi)->index();
      if (i >= _natoms)
        throw(TrajectoryError("updating group coords", _filename, "Atom index into trajectory frame is out of bounds"));
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
  }

  void AmberRst::seekFrameImpl(const uint i) {
  }

}
