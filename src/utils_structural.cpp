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
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>


#include <utils_structural.hpp>
#include <Coord.hpp>
#include <pdb.hpp>
#include <dcd.hpp>
#include <pdb_remarks.hpp>
#include <sfactories.hpp>

namespace loos {

  // Note: This *should* not throw a range_error, although in principal,
  //       the Remarks::operator[] does a range check and could...

  GCoord boxFromRemarks(const Remarks& r) {
    int n = r.size();
    int i;

    GCoord c(99999.99, 99999.99, 99999.99);

    for (i=0; i<n; i++) {
      std::string s = r[i];
      if (s.substr(0, 6) == " XTAL ") {
        std::stringstream is(s.substr(5));
        if (!(is >> c.x()))
          throw(ParseError("Unable to parse " + s));
        if (!(is >> c.y()))
          throw(ParseError("Unable to parse " + s));
        if (!(is >> c.z()))
          throw(ParseError("Unable to parse " + s));

        break;
      }
    }

    return(c);
  }



  // Note: This *should* not throw, although in principal,
  //       the Remarks::operator[] does a range check and could...

  bool remarksHasBox(const Remarks& r) {
    int n = r.size();
    for (int i = 0; i<n; i++) {
      std::string s = r[i];
      if (s.size() < 6)
	continue;
      if (s.substr(0, 6) == " XTAL ")
        return(true);
    }
    return(false);
  }



  AtomicGroup loadStructureWithCoords(const std::string& model_name, const std::string& coord_name) {
    AtomicGroup model = createSystem(model_name);
    if (!coord_name.empty()) {
      AtomicGroup coords = createSystem(coord_name);
      model.copyCoordinatesFrom(coords);
    }
    
    if (! model.hasCoords())
      throw(LOOSError("Error- no coordinates found in specified model(s)"));
    
    return(model);
  }
  
  AtomicGroup loadStructureWithCoords(const std::string& model_name, const std::string& type, const std::string& coord_name) {
    AtomicGroup model = createSystem(model_name, type);
    if (!coord_name.empty()) {
      AtomicGroup coords = createSystem(coord_name);
      model.copyCoordinatesFrom(coords);
    }
    
    if (! model.hasCoords())
      throw(LOOSError("Error- no coordinates found in specified model(s)"));
    
    return(model);
  }
  
  
  
  std::vector<uint> assignTrajectoryFrames(const pTraj& traj, const std::string& frame_index_spec, uint skip, uint stride)  {
    std::vector<uint> frames;
    
    if (frame_index_spec.empty())
      for (uint i=skip; i<traj->nframes(); i += stride)
        frames.push_back(i);
    else
      frames = parseRangeList<uint>(frame_index_spec, traj->nframes()-1);
    
    return(frames);
  }



};
