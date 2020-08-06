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






#if !defined(LOOS_HPP)
#define LOOS_HPP

// These are common system includes that nearly anybody who uses LOOS
// will probably be including anyway...

#include <iostream>
#include <iomanip>
#include <ios>
#include <sstream>
#include <fstream>

#include <ctime>
#include <cmath>

#include <string>
#include <vector>
#include <algorithm>

#include <stdexcept>

#include <cassert>

// These are the LOOS-specific includes...

#include <loos_defs.hpp>
#include <exceptions.hpp>
#include <utils.hpp>
#include <utils_random.hpp>
#include <utils_structural.hpp>

#include <Kernel.hpp>
#include <Parser.hpp>
#include <Selectors.hpp>


#include <Matrix44.hpp>
#include <XForm.hpp>
#include <Matrix.hpp>

#include <AtomicNumberDeducer.hpp>
#include <Atom.hpp>
#include <AtomicGroup.hpp>
#include <pdb.hpp>
#include <psf.hpp>
#include <amber.hpp>
#include <tinkerxyz.hpp>

#include <Trajectory.hpp>
#include <dcd.hpp>
#include <dcd_utils.hpp>
#include <MultiTraj.hpp>

#include <trajwriter.hpp>
#include <dcdwriter.hpp>
#include <xtcwriter.hpp>

#include <amber_traj.hpp>

#if defined(HAS_NETCDF)
#include <amber_netcdf.hpp>
#endif

#include <amber_rst.hpp>
#include <ccpdb.hpp>
#include <pdbtraj.hpp>
#include <tinker_arc.hpp>
#include <xtc.hpp>
#include <gro.hpp>
#include <trr.hpp>



#include <Geometry.hpp>
#include <ensembles.hpp>
#include <TimeSeries.hpp>

#include <Fmt.hpp>

#include <sfactories.hpp>

#include <loos_timer.hpp>
#include <ProgressCounters.hpp>
#include <ProgressTriggers.hpp>

#include <sorting.hpp>

#include <OptionsFramework.hpp>

#include <alignment.hpp>
#include <RnaSuite.hpp>
#endif


