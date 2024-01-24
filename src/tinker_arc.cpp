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

#include <tinker_arc.hpp>
#include <AtomicGroup.hpp>

namespace loos {

  void TinkerArc::init(void) {
    char buf[512];

    // Read the first frame to get the # of atoms...
    frame.read(*(ifs));
    _natoms = frame.size();
    indices.push_back(0l);
    cached_first = true;

    // Now determine the # of frames...
    while (1) {
      indices.push_back(ifs->tellg());
      ifs->getline(buf, sizeof(buf));
      if (frame.isPeriodic())
        ifs->getline(buf, sizeof(buf));
      for (uint i=0; i<_natoms; ++i)
        ifs->getline(buf, sizeof(buf));
      if (ifs->eof()) {
        break;
      }
    }

    _nframes = indices.size()-1;

    ifs->clear();
    ifs->seekg(indices[1]);
  }


  void TinkerArc::seekNextFrameImpl(void) {
    if (++current_index >= _nframes)
      at_end = true;
  }


  void TinkerArc::seekFrameImpl(const uint i) {
    if (i >= _nframes)
      throw(FileError(_filename, "Requested trajectory frame is out of range"));

    ifs->clear();
    ifs->seekg(indices[i]);
    if (ifs->fail())
      throw(FileError(_filename, "Cannot seek to the requested frame"));

    current_index = i;
    at_end = false;
  }


  bool TinkerArc::parseFrame(void) {
    if (ifs->eof() || at_end)
      return(false);

    // We embed the attempted read in a try catch
    // because the TinkerXYZ read() throws an error
    // if it can't read, so that the frame.size() test won't
    // keep us from walking off the end
    try {
      TinkerXYZ newframe;
      newframe.read(*(ifs));
      frame = newframe;
      if (frame.size() == 0) {
        at_end = true;
        return(false);
      }
    }
    catch(LOOSError& e) {
      return(false);
    }
    return(true);
  }


  std::vector<GCoord> TinkerArc::coords(void) const {
    std::vector<GCoord> result(_natoms);

    for (uint i=0; i<_natoms; i++)
      result[i] = frame[i]->coords();

    return(result);
  }


  void TinkerArc::updateGroupCoordsImpl(AtomicGroup& g)
  {
    for (AtomicGroup::iterator i = g.begin(); i != g.end(); ++i) {
      uint idx = (*i)->index();
      if (idx >= _natoms)
        throw(TrajectoryError("updating group coords", _filename, "Atom index into trajectory frame is out of bounds"));
      (*i)->coords(frame[idx]->coords());
    }

    // Handle periodic boundary conditions (if present)
    if (hasPeriodicBox()) {
      g.periodicBox(periodicBox());
    }
  }

}
