/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo, Alan Grossfield
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



#include <gro.hpp>



namespace loos {
  
  void Gromacs::read(std::istream& ifs) {
    std::string buf;

    getline(ifs, title_);

    // Get the # of atoms;;;
    getline(ifs, buf);
    int natoms = parseStringAs<int>(buf);
    while (natoms-- > 0) {
      getline(ifs, buf);
      
      int resid = parseStringAs<int>(buf, 0, 5);
      std::string resname = parseStringAs<std::string>(buf, 5, 5);
      std::string name = parseStringAs<std::string>(buf, 10, 5);
      int atomid = parseStringAs<int>(buf, 15, 5);
      float x = parseStringAs<float>(buf, 20, 8) * 10.0;
      float y = parseStringAs<float>(buf, 28, 8) * 10.0;
      float z = parseStringAs<float>(buf, 36, 8) * 10.0;
      
      // We ignore velocities...
      pAtom pa(new Atom);
      pa->index(_max_index++);
      pa->resid(resid);
      pa->id(atomid);
      pa->resname(resname);
      pa->name(name);
      pa->coords(GCoord(x,y,z));

      append(pa);
    }

    // Now process box...
    getline(ifs, buf);
    std::istringstream iss(buf);
    GCoord box;
    if (!(iss >> box[0] >> box[1] >> box[2]))
      throw(std::runtime_error("Cannot parse box '" + buf + "'"));
    periodicBox(box * 10.0);

    // Since the atomic field in .gro files is only 5-chars wide, it can
    // overflow.  if there are enough atoms to cause an overflow, manually
    // renumber everything...
    if (atoms.size() >= 100000)
      renumber();
  }

}
