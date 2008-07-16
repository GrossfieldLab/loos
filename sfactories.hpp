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






#if !defined(SFACTORIES_HPP)
#define SFACTORIES_HPP



#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <stdexcept>

#include <boost/algorithm/string.hpp>

#include <loos_defs.hpp>

#include <AtomicGroup.hpp>
#include <pdb.hpp>
#include <psf.hpp>
#include <amber.hpp>

#include <Trajectory.hpp>
#include <dcd.hpp>
#include <amber_traj.hpp>

using namespace std;
using namespace boost;


namespace loos {
  AtomicGroup createSystem(const string&);
  pTraj createTrajectory(const string&, const AtomicGroup&);
};


#endif
