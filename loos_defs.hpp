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






#if !defined(LOOSDEFS_HPP)
#define LOOSDEFS_HPP


#include <Coord.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

namespace loos {

  typedef double greal;
  typedef long gint;

  typedef float dcd_real;
  typedef double dcd_double;

  typedef Coord<double> GCoord;
  typedef boost::shared_ptr<GCoord> pGCoord;

  // Trajectory and subclasses...
  class Trajectory;
  class DCD;
  class AmberTraj;
  class CCPDB;
  class TinkerArc;
  class PDBTraj;


  typedef boost::shared_ptr<Trajectory> pTraj;
  typedef boost::shared_ptr<DCD> pDCD;
  typedef boost::shared_ptr<AmberTraj> pAmberTraj;
  typedef boost::shared_ptr<CCPDB> pCCPDB;
  typedef boost::shared_ptr<TinkerArc> pTinkerArc;
  typedef boost::shared_ptr<PDBTraj> pPDBTraj;

  // AtomicGroup and subclasses (i.e. systems formats)
  class AtomicGroup;
  class PDB;
  class PSF;
  class Amber;
  class TinkerXYZ;

  typedef boost::shared_ptr<AtomicGroup> pAtomicGroup;
  typedef boost::shared_ptr<PDB> pPDB;
  typedef boost::shared_ptr<PSF> pPSF;
  typedef boost::shared_ptr<Amber> pAmber;
  typedef boost::shared_ptr<TinkerXYZ> pTinkerXYZ;


  const uint kilobytes = 1024;
  const uint megabytes = kilobytes * kilobytes;
  const uint gigabytes = megabytes * kilobytes;

}

#endif


