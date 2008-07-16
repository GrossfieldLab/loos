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


#include <sfactories.hpp>

AtomicGroup loos::createSystem(const string& s) {

  if (iends_with(s, ".pdb")) {
    PDB pdb(s);
    return(pdb);
  } else if (iends_with(s, ".psf")) {
    PSF psf(s);
    return(psf);
  } else if (iends_with(s, ".prmtop")) {
    Amber amber(s);
    return(amber);
  } else
    throw(runtime_error("Error- cannot divine file type from name '" + s + "'"));
}



pTraj loos::createTrajectory(const string& s, const AtomicGroup& g) {

  if (iends_with(s, ".dcd")) {
    pDCD pd(new DCD(s));
    pTraj pt(pd);
    return(pt);
  } else if (iends_with(s, ".mdcrd")) {
    pAmberTraj pat(new AmberTraj(s, g.size()));
    pTraj pt(pat);
    return(pt);
  } else
    throw(runtime_error("Error- cannot divine file type from name '" + s + "'"));
}

