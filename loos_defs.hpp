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

#if __GNUC__ < 4
#error LOOS Requires GCC-4.0.1 or higher
#endif

#include <sys/types.h>


#if defined(REQUIRES_UINT)
typedef unsigned int     uint;
#endif

#if defined(REQUIRES_ULONG)
typedef unsigned long    ulong;
#endif

#if defined(__linux__)

extern "C" {

#include <cblas.h>
#include <clapack.h>

  void dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);
  void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
  void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
  void dggev_(char*, char*, int*, double*, int*, double*, int*, double*, double*, double*, double*, int*, double*, int*, double*, int*, int*);
  void dpotri_(char*, int*, double*, int*, int*);

  void sgesvd_(char*, char*, int*, int*, float*, int*, float*, float*, int*, float*, int*, float*, int*, int*);
  void sgemm_(char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
  void sggev_(char*, char*, int*, float*, int*, float*, int*, float*, float*, float*, float*, int*, float*, int*, float*, int*, int*);
  void ssyev_(char*, char*, int*, float*, int*, float*, float*, int*, int*);
  void ssygv_(int*, char*, char*, int*, float*, int*, float*, int*, float*, float*, int*, int*);
  void ssygvx_(int*, char*, char*, char*, int*, float*, int*, float*, int*, float*, float*, int*, int*, float*, int*, float*, float*, int*, float*, int*, int*, int*, int*);
  double dlamch_(char*);
}


typedef int f77int;

#elif defined(__APPLE__)

#include <vecLib/vecLib.h>

typedef __CLPK_integer f77int;

#else

#error You are building in an unsupported environment.  You will need to specify how to access ATLAS.

#endif




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
  class Atom;
  class Trajectory;
  class DCD;
  class AmberTraj;
  class CCPDB;
  class TinkerArc;
  class PDBTraj;
  class XTC;
  class TRR;


  typedef boost::shared_ptr<Atom> pAtom;
  typedef boost::shared_ptr<Trajectory> pTraj;
  typedef boost::shared_ptr<DCD> pDCD;
  typedef boost::shared_ptr<AmberTraj> pAmberTraj;
  typedef boost::shared_ptr<CCPDB> pCCPDB;
  typedef boost::shared_ptr<TinkerArc> pTinkerArc;
  typedef boost::shared_ptr<PDBTraj> pPDBTraj;
  typedef boost::shared_ptr<XTC> pXTC;
  typedef boost::shared_ptr<TRR> pTRR;

  // AtomicGroup and subclasses (i.e. systems formats)
  class AtomicGroup;
  class PDB;
  class PSF;
  class Amber;
  class AmberRst;
  class TinkerXYZ;
  class Gromacs;


  typedef boost::shared_ptr<AtomicGroup> pAtomicGroup;
  typedef boost::shared_ptr<PDB> pPDB;
  typedef boost::shared_ptr<PSF> pPSF;
  typedef boost::shared_ptr<Amber> pAmber;
  typedef boost::shared_ptr<AmberRst> pAmberRst;
  typedef boost::shared_ptr<TinkerXYZ> pTinkerXYZ;
  typedef boost::shared_ptr<Gromacs> pGromacs;


  const uint kilobytes = 1024;
  const uint megabytes = kilobytes * kilobytes;
  const uint gigabytes = megabytes * kilobytes;

}

#endif


