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





%header %{
#include <iostream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <exception>
#include <stdexcept>


#include <boost/random.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include <ctime>


#include <loos_defs.hpp>
#include <exceptions.hpp>
#include <Coord.hpp>
#include <pdb_remarks.hpp>

#include <utils.hpp>
%}


//! Namespace for most things not already encapsulated within a class.
namespace loos {

  std::string findBaseName(const std::string&);
  std::string getNextLine(std::istream& is, int* lineno);
  AtomicGroup selectAtoms(const AtomicGroup&, const std::string) throw (loos::NullResult);
  AtomicGroup loadStructureWithCoords(const std::string& model, const std::string& cooords);
  std::vector<uint> assignTrajectoryFrames(const pTraj& traj, const std::string& frame_index_spec, uint skip = 0, uint stride = 1);
};
