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

#include <psf.hpp>
#include <exceptions.hpp>


namespace loos {


  PSF* PSF::clone(void) const {
    return(new PSF(*this));
  }

  PSF PSF::copy(void) const {
    AtomicGroup grp = this->AtomicGroup::copy();
    PSF p(grp);
    
    // Add PSF specific member data copies here...
    return(p);
  }



  void PSF::read(std::istream& is) {
    std::string input;

    // first line is the PSF header
    if (!getline(is, input))
      throw(std::runtime_error("Failed reading first line of psf"));
    if (input.substr(0,3) != "PSF")
      throw(std::runtime_error("PSF detected a non-PSF file"));

    // second line is blank
    if (!getline(is, input))
      throw(std::runtime_error("PSF failed reading first header blank"));

    // third line is title header
    getline(is, input);
    int num_title_lines;
    if (!(std::stringstream(input) >> num_title_lines)) 
      throw(std::runtime_error("PSF has malformed title header"));
        
    // skip the rest of the title
    for (int i=0; i<num_title_lines; i++) 
      getline(is, input);
        
    // verify nothing went wrong
    if (!(is.good()))
      // Yes, I know, I should figure out what went wrong instead
      // of running home crying.  Sorry, Tod...
      throw(std::runtime_error("PSF choked reading the header"));

    // next line is blank 
    if (!getline(is, input))
      throw(std::runtime_error("PSF failed reading second header blank"));

    // next line is the number of atoms
    
    if (!getline(is, input))
      throw(std::runtime_error("PSF failed reading natom line"));
    int num_atoms;
    if (!(std::stringstream(input) >> num_atoms))
      throw(std::runtime_error("PSF has malformed natom line"));

    for (int i=0; i<num_atoms; i++) {
      if (!getline(is, input))
        throw(std::runtime_error("Failed reading PSF atom line "));
      //throw(std::runtime_error("Failed reading PSF atom line " + std::string(i)));
      parseAtomRecord(input);
    }

    // next line is blank 
    if (!getline(is, input))
      throw(std::runtime_error("PSF failed reading blank after atom lines"));

    // next block of lines is the list of bonds
    // Bond title line
    if (!getline(is, input))
      throw(std::runtime_error("PSF failed reading nbond line"));
    int num_bonds;
    if (!(std::stringstream(input) >> num_bonds))
      throw(std::runtime_error("PSF has malformed nbond line"));

    int bonds_found = 0;
    getline(is, input);
    while (input.size() > 0) { // end of the block is marked by a blank line
      int ind1, ind2;
      std::stringstream s(input);
      std::string sind1, sind2;   // Hold string representations of id's before we hybrid36 decode them...
      while (s.good()) {
        s >> sind1;
        s >> sind2;

        ind1 = parseStringAsHybrid36(sind1);
        ind2 = parseStringAsHybrid36(sind2);

        if (ind1 > num_atoms || ind2 > num_atoms)
          throw(LOOSError("PSF bond error: bound atomid exceeds number of atoms."));

        ind1--;  // our indices are 1 off from the numbering in the pdb/psf file
        ind2--;
        pAtom pa1 = getAtom(ind1);                
        pAtom pa2 = getAtom(ind2);                
        pa1->addBond(pa2);
        pa2->addBond(pa1);
        bonds_found++; 
      }
      getline(is, input);
    }
    // sanity check
    if (bonds_found != num_bonds) 
      throw(std::runtime_error("PSF number of bonds disagrees with number found"));

    deduceAtomicNumberFromMass();
    setGroupConnectivity();
  }



  void PSF::parseAtomRecord(const std::string s) {
    gint index;
    std::string segname;
    gint resid;
    std::string resname;
    std::string atomname;
    std::string atomtype;
    greal charge;
    greal mass;
    gint fixed;

    pAtom pa(new Atom);
         
    std::stringstream ss(s);

    // A hack to handle hybrid-36 resids & atomids
    std::string buf;
    ss >> buf;
    index = parseStringAsHybrid36(buf);
    pa->id(index);
     
    ss >> segname;
    pa->segid(segname);


    ss >> buf;
    resid = parseStringAsHybrid36(buf);
    pa->resid(resid);

    ss >> resname;
    pa->resname(resname);

    ss >> atomname;
    pa->name(atomname);

    // If this is a charmm psf, the atomtype will be an integer.
    // NAMD/XPLOR psfs use the symbolic atomtype, which must start with a letter
    // At the moment, the Atom class doesn't care about this value (it's mostly
    // used in charmm and namd as a means to look up parameters), so we're going to
    // discard it.  However, if we ever decide we're going to use this, we'll need
    // to keep track of the distinction between charmm and namd usage.
    ss >> atomtype;

    ss >> charge;
    pa->charge(charge);

    ss >> mass;
    pa->mass(mass);

    // Is the atom fixed or mobile?
    // for now, we're going to silently drop this
    ss >> fixed;  

    append(pa);
  }

}
