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
#include <sys/stat.h>


#include <boost/algorithm/string.hpp>

#include <AtomicGroup.hpp>
#include <pdb.hpp>
#include <psf.hpp>
#include <amber.hpp>

#include <Trajectory.hpp>
#include <dcd.hpp>
#include <amber_traj.hpp>
#include <amber_rst.hpp>
#include <ccpdb.hpp>
#include <charmm.hpp>
#include <tinkerxyz.hpp>
#include <tinker_arc.hpp>
#include <gro.hpp>
#include <xtc.hpp>
#include <trr.hpp>


namespace loos {

  std::string availableSystemFileTypes() {
    std::string types =
      "psf\tCHARMm/NAMD PSF\n"
      "pdb\tCHARMm/NAMD PDB\n"
      "crd\tCHARMm CRD\n"
      "xyz\tTinker XYZ\n"
      "prtmop\tAmber PRMTOP\n"
      "gro\tGROMACS GRO\n";
    return(types);
  }


  pAtomicGroup createSystemPtr(const std::string& filename, const std::string& filetype) {

    pAtomicGroup pag;

    if (filetype == "pdb") {
      pPDB p(new PDB(filename));
      pag = p;
    } else if (filetype == "psf") {
      pPSF p(new PSF(filename));
      pag = p;
    } else if (filetype == "prmtop") { 
      pAmber p(new Amber(filename));
      pag = p;
    } else if (filetype == "xyz") {
      pTinkerXYZ p(new TinkerXYZ(filename));
      pag = p;
    } else if (filetype == "gro") {
      pGromacs p(new Gromacs(filename));
      pag = p;

    } else
      throw(std::runtime_error("Error- unknown system file type '" + filetype + "' for file '" + filename + "'"));

    return(pag);
  }



  pAtomicGroup createSystemPtr(const std::string& filename) {

    size_t extension_pos = filename.rfind('.');
    if (extension_pos == filename.npos)
      throw(std::runtime_error("Error- system filename must end in an extension or the filetype must be explicitly specified"));

    std::string filetype = filename.substr(extension_pos+1);
    boost::to_lower(filetype);
    return(createSystemPtr(filename, filetype));
  }


  AtomicGroup createSystem(const std::string& filename) {
    return(*(createSystemPtr(filename)));
  }

  AtomicGroup createSystem(const std::string& filename, const std::string& filetype) {
    return(*(createSystemPtr(filename, filetype)));
  }



  std::string availableTrajectoryFileTypes() {
    std::string types =
      "dcd\tCHARMm/NAMD DCD\n"
      "mdcrd\tAmber MDCRD\n"
      "rst\tAmber RST\n"
      "rst7\tAmber RST\n"
      "inpcrd\tAmber RST\n"
      "pdb\tConcatenated CHARMm/NAMD PDB Trajectory\n"
      "arc\tTinker ARC\n"
      "xtc\tGROMACS XTC\n"
      "trr\tGROMACS TRR\n";
    return(types);
  }


  pTraj createTrajectory(const std::string& filename, const std::string& filetype, const AtomicGroup& g) {
    
    if (filetype == "dcd") {
      pDCD pd(new DCD(filename));
      pTraj pt(pd);
      return(pt);
    } else if (filetype == "mdcrd") {
      pAmberTraj pat(new AmberTraj(filename, g.size()));
      pTraj pt(pat);
      return(pt);
    } else if (filetype == "rst"
               || filetype == "rst7"
               || filetype == "inpcrd") {
      pAmberRst par(new AmberRst(filename, g.size()));
      pTraj pt(par);
      return(pt);
    } else if (filetype == "pdb") {
      pCCPDB ppdb(new CCPDB(filename));
      pTraj pt(ppdb);
      return(pt);
    } else if (filetype == "arc") {
      pTinkerArc pta(new TinkerArc(filename));
      pTraj pt(pta);
      return(pt);
    } else if (filetype == "xtc") {
      pXTC pxtc(new XTC(filename));
      pTraj pt(pxtc);
      return(pt);
    } else if (filetype == "trr") {
      pTRR ptrr(new TRR(filename));
      pTraj pt(ptrr);
      return(pt);

    } else
      throw(std::runtime_error("Error- unknown trajectory file type '" + filetype + "' for file '" + filename + "'"));
  }


  pTraj createTrajectory(const std::string& filename, const AtomicGroup& g) {
    size_t extension_pos = filename.rfind('.');
    if (extension_pos == filename.npos)
      throw(std::runtime_error("Error- trajectory filename must end in an extension or the filetype must be explicitly specified"));

    std::string filetype = filename.substr(extension_pos+1);
    boost::to_lower(filetype);
    return(createTrajectory(filename, filetype, g));
  }

}



