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

#include <tinkerxyz.hpp>
#include <Atom.hpp>

namespace loos {


  TinkerXYZ* TinkerXYZ::clone(void) const {
    return(new TinkerXYZ(*this));
  }

  TinkerXYZ TinkerXYZ::copy(void) const {
    AtomicGroup grp = this->AtomicGroup::copy();
    TinkerXYZ p(grp);
    
    // Add TinkerXYZ specific member data copies here...
    return(p);
  }


  void TinkerXYZ::read(std::istream& is) {
    std::string input;

    // first line is the header, first field is number of atoms
    if (!getline(is, input))
      throw(FileReadError(_filename, "Failed reading first line of xyz"));
    int num_atoms = 0;
    if (!(std::stringstream(input) >> num_atoms))
      throw(FileReadError(_filename, "TinkerXYZ has malformed header"));

    if (!getline(is, input)) 
      throw(FileReadError(_filename, "Failed reading 2nd line of Tinker file"));
    if (parseBoxRecord(input)) {
      if (!getline(is, input)) {
        throw(FileReadError(_filename, "Failed reading first atom line of TinkerXYZ"));
      }
    }



    // Read the lines
    for (int i=0; i<num_atoms; i++) {
      parseAtomRecord(input);
      if ((i < (num_atoms-1)) && (!getline(is, input)))
        throw(FileReadError(_filename, "Failed reading TinkerXYZ atom line "));
    }

    // Assume all XYZ files have connectivity
    setGroupConnectivity();
  }


  bool TinkerXYZ::parseBoxRecord(const std::string& s) {
     greal x,y,z;
     std::string tmp;
     std::stringstream ss(s);
     ss >> x;
     ss >> tmp;
     if (isalpha(tmp[0])) {
       // this is an atom record, not a box record
       return false;
     }
     else {
       y = atof(tmp.c_str());
     }
     ss >> z;

     periodicBox(x,y,z);
     return true;
  }


  void TinkerXYZ::parseAtomRecord(const std::string& s) {


    gint index;
    //std::string segname("");       // Tinker doesn't have segments
    //gint resid=1;             // Tinker doesn't have residues
    //std::string resname("");       // Tinker doesn't have residues
    std::string atomname;  
    //greal charge=0.0;
    //greal mass=1.0;
    //gint atomic_number = 1;

    pAtom pa(new Atom);
    pa->index(_max_index++);
         
    std::stringstream ss(s);

    ss >> index;
    pa->id(index);
     
    ss >> atomname;
    pa->name(atomname);

    greal x,y,z;
    ss >> x;
    ss >> y;
    ss >> z;
    pa->coords(GCoord(x,y,z));

    int atomtype;
    ss >> atomtype;
    pa->atomType(atomtype);

    // Now read in the atoms to which this atom is bonded
    int bonded_atom;
    while (ss >> bonded_atom)
      {
        // Probably should verify this value is sane
        pa->addBond(bonded_atom); 
      }

    // Tinker XYZ files don't have segments or residues, so these
    // properties get set to the class defaults by the constructor
    //pa->segid(segname);  
    //pa->resid(resid);
    //pa->resname(resname);

    // TODO: once we support reading Tinker param files, we should set these
    // to their correct values
    //pa->charge(charge);
    //pa->mass(mass);
    //pa->atomic_number(atomic_number);


    append(pa);
  }

}
