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




%module sfactories
%include "std_string.i"
%include "AtomicGroup.i"
%include "Trajectory.i"

%{
#include <loos_defs.hpp>
  //  typedef boost::shared_ptr<Trajectory>   pTraj;
%}

namespace loos {
  //! Factory function for reading in structure files.
  /*!
   * This function will try to determine the filetype for a structure
   * file by examining the suffix of the file.  It will return an
   * AtomicGroup copy of the input structure.
   *
   * Amber files are handled differently...  If the filename ends with
   * ".prmtop", then LOOS will check for a corresponding ".inpcrd"
   * file.  If one is found, then the coords will be loaded into the
   * group.  Otherwise, the prmtop will be loaded without coords and
   * returned.
   *
   */
  AtomicGroup createSystem(const std::string&);

  //pAtomicGroup createSystemPtr(const std::string&);

  //! Factory function for reading in a trajectory file.
  /*!
   * This function will try to determine the filetype for a trajectory
   * by examining the suffix of the file.  It will return a boost
   * shared pointer to a new Trajectory object.
   *
   * It is \e very \e important to understand that the object returned
   * by this function must behave polymorphically.  That's why it is
   * wrapped in a boost shared pointer to the base class.  Do not
   * try to deference it and assign it to a Trajectory object...
   */
  pTraj createTrajectory(const std::string&, const AtomicGroup&);
};

