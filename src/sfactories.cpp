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


#include <string>

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

#include <amber_netcdf.hpp>
#include <amber_rst.hpp>
#include <ccpdb.hpp>
#include <charmm.hpp>
#include <tinkerxyz.hpp>
#include <tinker_arc.hpp>
#include <gro.hpp>
#include <xtc.hpp>
#include <trr.hpp>


#include <trajwriter.hpp>
#include <dcdwriter.hpp>
#include <xtcwriter.hpp>

namespace loos {


  namespace internal {
    struct SystemNameBindingType {
      std::string suffix;
      std::string type;
      pAtomicGroup (*creator)(const std::string& fname);
    };

    SystemNameBindingType system_name_bindings[] = {
      { "prmtop", "Amber", &Amber::create },
      { "crd", "CHARMM CRD", &CHARMM::create },
      { "pdb", "CHARMM/NAMD PDB", &PDB::create },
      { "psf", "CHARMM/NAMD PSF", &PSF::create },
      { "gro", "Gromacs", &Gromacs::create },
      { "xyz", "Tinker", &TinkerXYZ::create },
      { "", "", 0}
    };
  }




  std::string availableSystemFileTypes(const std::string& prefix) {
    std::string types;
    for (internal::SystemNameBindingType* p = internal::system_name_bindings; p->creator != 0; ++p) {
      types += prefix + p->suffix + " = " + p->type + "\n";
    }

    return(types);
  }


  pAtomicGroup createSystemPtr(const std::string& filename, const std::string& filetype) {

    for (internal::SystemNameBindingType* p = internal::system_name_bindings; p->creator != 0; ++p)
      if (p->suffix == filetype)
        return(*(p->creator))(filename);

    throw(std::runtime_error("Error- unknown system file type '" + filetype + "' for file '" + filename + "'.  Try --help to see available types."));
  }



  pAtomicGroup createSystemPtr(const std::string& filename) {

    boost::tuple<std::string, std::string> names = splitFilename(filename);
    std::string suffix = boost::get<1>(names);
    if (suffix.empty())
      throw(std::runtime_error("Error- system filename must end in an extension or the filetype must be explicitly specified"));

    boost::to_lower(suffix);
    return(createSystemPtr(filename, suffix));
  }


  AtomicGroup createSystem(const std::string& filename) {
    return(*(createSystemPtr(filename)));
  }

  AtomicGroup createSystem(const std::string& filename, const std::string& filetype) {
    return(*(createSystemPtr(filename, filetype)));
  }


  namespace internal {
    struct TrajectoryNameBindingType {
      std::string suffix;
      std::string type;
      pTraj (*creator)(const std::string& fname, const AtomicGroup& model);
    };

    TrajectoryNameBindingType trajectory_name_bindings[] = {
      { "crd", "Amber Traj (NetCDF/Amber)", &AmberNetcdf::create},
      { "mdcrd", "Amber Traj (NetCDF/Amber)", &AmberNetcdf::create},
      { "nc", "Amber Traj (NetCDF)", &AmberNetcdf::create},
      { "netcdf", "Amber Traj (NetCDF)", &AmberNetcdf::create},
      { "inpcrd", "Amber Restart", &AmberRst::create},
      { "rst", "Amber Restart", &AmberRst::create},
      { "rst7", "Amber Restart", &AmberRst::create},
      { "dcd", "CHARMM/NAMD DCD", &DCD::create},
      { "pdb", "Concatenated PDB", &CCPDB::create},
      { "trr", "Gromacs TRR", &TRR::create},
      { "xtc", "Gromacs XTC", &XTC::create},
      { "arc", "Tinker ARC", &TinkerArc::create},
      { "", "", 0}
    };



  }


  std::string availableTrajectoryFileTypes(const std::string& prefix) {
    std::string types;
    for (internal::TrajectoryNameBindingType* p = internal::trajectory_name_bindings; p->creator != 0; ++p) {
      types += prefix + p->suffix + " = " + p->type + "\n";
    }

    return(types);

  }


  pTraj createTrajectory(const std::string& filename, const std::string& filetype, const AtomicGroup& g) {

    // First, check to make sure AtomicGroup has index information...
    if (!g.allHaveProperty(Atom::indexbit))
      throw(LOOSError("Model passed to createTrajectory() does not have atom index information."));

    for (internal::TrajectoryNameBindingType* p = internal::trajectory_name_bindings; p->creator != 0; ++p) {
      if (p->suffix == filetype) {
        return((*(p->creator))(filename, g) );
      }
    }
    throw(std::runtime_error("Error- unknown input trajectory file type '" + filetype + "' for file '" + filename + "'.  Try --help to see available types."));

  }


  pTraj createTrajectory(const std::string& filename, const AtomicGroup& g) {
    boost::tuple<std::string, std::string> names = splitFilename(filename);
    std::string suffix = boost::get<1>(names);

    if (suffix.empty())
      throw(std::runtime_error("Error- trajectory filename must end in an extension or the filetype must be explicitly specified"));

    boost::to_lower(suffix);
    return(createTrajectory(filename, suffix, g));
  }


  namespace internal {
    struct OutputTrajectoryNameBindingType {
      std::string suffix;
      std::string type;
      pTrajectoryWriter (*creator)(const std::string& fname, const bool append);
    };

    OutputTrajectoryNameBindingType output_trajectory_name_bindings[] = {
      { "dcd", "NAMD DCD", &DCDWriter::create},
      { "xtc", "Gromacs XTC (compressed trajectory)", &XTCWriter::create},
      { "", "", 0}
    };

  }


  std::string availableOutputTrajectoryFileTypes(const std::string& prefix) {
    std::string types;
    for (internal::OutputTrajectoryNameBindingType* p = internal::output_trajectory_name_bindings; p->creator != 0; ++p) {
      types += prefix + p->suffix + " = " + p->type + "\n";
    }

    return(types);
  }


  pTrajectoryWriter createOutputTrajectory(const std::string& filename, const std::string& type, const bool append) {

    for (internal::OutputTrajectoryNameBindingType* p = internal::output_trajectory_name_bindings; p->creator != 0; ++p) {
      if (p->suffix == type) {
        return((*(p->creator))(filename, append));
      }
    }

    throw(std::runtime_error("Error- unknown output trajectory file type '" + type + "' for file '" + filename + "'.  Try --help to see available types."));
  }


  pTrajectoryWriter createOutputTrajectory(const std::string& filename, const bool append) {
    boost::tuple<std::string, std::string> names = splitFilename(filename);
    std::string suffix = boost::get<1>(names);
    if (suffix.empty())
      throw(std::runtime_error("Error- output trajectory filename must end in an extension or the filetype must be explicitly specified"));

    return(createOutputTrajectory(filename, suffix, append));
  }


}
