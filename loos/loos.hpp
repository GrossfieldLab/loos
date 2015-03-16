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

#include <loos/loos_defs.hpp>
#include <loos/exceptions.hpp>
#include <loos/utils.hpp>
#include <loos/utils_random.hpp>
#include <loos/utils_structural.hpp>

#include <loos/Kernel.hpp>
#include <loos/Parser.hpp>
#include <loos/Selectors.hpp>


#include <loos/Matrix44.hpp>
#include <loos/XForm.hpp>
#include <loos/Matrix.hpp>

#include <loos/AtomicNumberDeducer.hpp>
#include <loos/Atom.hpp>
#include <loos/AtomicGroup.hpp>
#include <loos/pdb.hpp>
#include <loos/psf.hpp>
#include <loos/amber.hpp>
#include <loos/tinkerxyz.hpp>

#include <loos/Trajectory.hpp>
#include <loos/dcd.hpp>
#include <loos/dcd_utils.hpp>

#include <loos/trajwriter.hpp>
#include <loos/dcdwriter.hpp>
#include <loos/xtcwriter.hpp>

#include <loos/amber_traj.hpp>

#if defined(HAS_NETCDF)
#include <loos/amber_netcdf.hpp>
#endif

#include <loos/amber_rst.hpp>
#include <loos/ccpdb.hpp>
#include <loos/pdbtraj.hpp>
#include <loos/tinker_arc.hpp>
#include <loos/xtc.hpp>
#include <loos/gro.hpp>
#include <loos/trr.loos/hpp>



#include <loos/Geometry.hpp>
#include <loos/ensembles.hpp>
#include <loos/TimeSeries.hpp>

#include <loos/Fmt.hpp>

#include <loos/sfactories.hpp>

#include <loos/timer.hpp>
#include <loos/ProgressCounters.hpp>
#include <loos/ProgressTriggers.hpp>

#include <loos/sorting.hpp>

#include <loos/OptionsFramework.hpp>

#endif


