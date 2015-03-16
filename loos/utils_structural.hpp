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




#if !defined(LOOS_UTILS_STRUCTURAL_HPP)
#define LOOS_UTILS_STRUCTURAL_HPP

#include <vector>

#include <loos/loos_defs.hpp>
#include <loos/exceptions.hpp>

namespace loos {
  //! Extract the Alan-style box-size from a PDB Remarks block.
  /** Returns a GCoord(99999.99, 99999.99, 99999.99) if there is no box
   *  info found in the remarks block.
   */
  GCoord boxFromRemarks(const Remarks&);

  //! Checks to see if a Remarks block has an Alan-style box size in it.
  bool remarksHasBox(const Remarks&);

  //! Loads a structure and optional coordinates
  AtomicGroup loadStructureWithCoords(const std::string& model, const std::string& cooords);

  AtomicGroup loadStructureWithCoords(const std::string& model, const std::string& type, const std::string& cooords);

  //! Builds a list of trajectory indices (frame_index_spec supercedes skip)
  std::vector<uint> assignTrajectoryFrames(const pTraj& traj, const std::string& frame_index_spec, uint skip = 0, uint stride = 1);

};

#endif


