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

#include <charmm.hpp>
#include <utils.hpp>
#include <boost/algorithm/string.hpp>




namespace loos {


  CHARMM* CHARMM::clone(void) const {
    return(new CHARMM(*this));
  }

  CHARMM CHARMM::copy(void) const {
    AtomicGroup grp = this->AtomicGroup::copy();
    CHARMM p(grp);
    
    return(p);
  }



  //! The file starts with one or more comment lines, which begin with a "*"
  //! These will be discarded (for now at least -- I suppose I could mimic
  //! what's done with a PDB file and store them, but I don't feel like it 
  //! now.
  //! Note: For now, we ignore the "RESID" field, but stuff the weighting
  //! field into the occupancy.
  void CHARMM::read(std::istream& is) {
    std::string input;


    getline(is, input);
    while (input[0] == '*') {
        getline(is, input);
        if (!is.good()) {
            throw(std::runtime_error("Failed reading CHARMM header"));
        }
    }

    // next line is just the number of atoms, by itself
    unsigned int num_atoms = atoi(input.c_str());

    // now loop and read the coordinates
    // Note: there are two different formats, depending on the number of atoms
    int atom_num = -1;
    int res_num = -1;
    std::string res_name;
    std::string atom_name;
    float x,y,z;
    std::string segid;
    float weight;
    
    for (unsigned int i=1; i<= num_atoms; ++i)  {
        getline(is, input);
        if (!is.good()) {
            throw(std::runtime_error("Failed reading CHARMM coordinates"));
        }
        if (num_atoms < 100000) {
            atom_num = parseStringAs<int>(input, 0, 5);
            res_num = parseStringAs<int>(input,5,5);
            res_name = parseStringAs<std::string>(input,11,4);
            atom_name = parseStringAs<std::string>(input,16,4);
            x = parseStringAs<float>(input,20,10);
            y = parseStringAs<float>(input,30,10);
            z = parseStringAs<float>(input,40,10);
            segid = parseStringAs<std::string>(input,51,4);
            weight = parseStringAs<float>(input,63,10);
        }
        else {
            atom_num = parseStringAs<int>(input, 0, 10);
            res_num = parseStringAs<int>(input,10,10);
            res_name = parseStringAs<std::string>(input,22,8);
            atom_name = parseStringAs<std::string>(input,32,8);
            x = parseStringAs<float>(input,40,20);
            y = parseStringAs<float>(input,60,20);
            z = parseStringAs<float>(input,80,20);
            segid = parseStringAs<std::string>(input,82,8);
            weight = parseStringAs<float>(input,100,20);
        }

        // Create a new atom and fill in the values
        pAtom pa(new Atom);
        pa->id(atom_num);
        pa->resid(res_num);
        pa->name(atom_name);
        pa->resname(res_name);
        pa->coords().x() = x;
        pa->coords().y() = y;
        pa->coords().z() = z;
        pa->occupancy(weight);

        // add the new atom to the CHARMM list
        append(pa);          
    }
  }
        

}
