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


#include <ios>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include <ctype.h>

#include "loos.hpp"
#include "psf.hpp"

void PSF::read(istream& is) {
    string input;

    // first line is the PSF header
    if (!getline(is, input))
        throw(runtime_error("Failed reading first line of psf"));
    if (input.substr(0,3) != "PSF")
        throw(runtime_error("PSF detected a non-PSF file"));

    // second line is blank
    if (!getline(is, input))
        throw(runtime_error("PSF failed reading first header blank"));

    // third line is title header
    getline(is, input);
    int num_title_lines;
    if (!(stringstream(input) >> num_title_lines)) 
        throw(runtime_error("PSF has malformed title header"));
        
    // skip the rest of the title
    for (int i=0; i<num_title_lines; i++) 
        getline(is, input);
        
    // verify nothing went wrong
    if (!(is.good()))
        // Yes, I know, I should figure out what went wrong instead
        // of running home crying.  Sorry, Tod...
        throw(runtime_error("PSF choked reading the header"));

    // next line is blank 
    if (!getline(is, input))
        throw(runtime_error("PSF failed reading second header blank"));

    // next line is the number of atoms
    
    if (!getline(is, input))
        throw(runtime_error("PSF failed reading natom line"));
    int num_atoms;
    if (!(stringstream(input) >> num_atoms))
        throw(runtime_error("PSF has malformed natom line"));

    for (int i=0; i<num_atoms; i++) {
        if (!getline(is, input))
            throw(runtime_error("Failed reading PSF atom line "));
            //throw(runtime_error("Failed reading PSF atom line " + string(i)));
        parseAtomRecord(input);
    }

    // next line is blank 
    if (!getline(is, input))
        throw(runtime_error("PSF failed reading blank after atom lines"));

    // next block of lines is the list of bonds
    // Bond title line
    if (!getline(is, input))
        throw(runtime_error("PSF failed reading nbond line"));
    int num_bonds;
    if (!(stringstream(input) >> num_bonds))
        throw(runtime_error("PSF has malformed nbond line"));

    int bonds_found = 0;
    getline(is, input);
    while (input.size() > 0) { // end of the block is marked by a blank line
        int ind1, ind2;
        stringstream s(input);
        while (s.good()) {
            s >> ind1;
            s >> ind2;
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
        throw(runtime_error("PSF number of bonds disagrees with number found"));

}



void PSF::parseAtomRecord(const string s) {
    gint index;
    string segname;
    gint resid;
    string resname;
    string atomname;
    string atomtype;
    greal charge;
    greal mass;
    gint fixed;

    pAtom pa(new Atom);
         
    stringstream ss(s);

    ss >> index;
    pa->id(index);
     
    ss >> segname;
    pa->segid(segname);

    ss >> resid;
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
    pa->atomic_number(deduceAtomicNumber(pa));

    // Is the atom fixed or mobile?
    // for now, we're going to silently drop this
    ss >> fixed;  

    // now intialize the rest of the stuff
    pa->recordName(string(""));

#if 0
    // Removed -- this is no longer necessary because we have the 
    // bitmask which lets us mark the coordinates as unset.  In fact
    // this way would be broken, because it would cause the coordinates 
    // bit to be set to true
    GCoord c= GCoord(99999.99, 99999.99, 99999.99); // Marks
						    // uninitialized
						    // coords 
    pa->coords(c);
#endif 

    append(pa);
}


/** This is a bit ad-hoc, and doesn't include all possible nuclei, but it
 *  covers everything that's likely to show up in a biomolecular simulation.
 *  However, at the moment it doesn't handle isotopes (eg deuterium).
 *  I got this list from the standard CHARMM topology files -- obviously, more
 *  atom types (eg other cations) could be added as needed.  The ranges I used
 *  are kind of large, since I recall the different force fields used to disagree
 *  about some of the masses.
 */
int PSF::deduceAtomicNumber(pAtom pa) {
    double mass = pa->mass();
    int an=-1;
    if (mass < 1.0)
        throw(out_of_range("Atomic mass less than 1.0 in psf"));
    else if (mass < 1.1) // Hydrogen = 1.0080
        an=1;
    else if ( (mass >= 4.0) && (mass <= 4.1) )  // Helium = 4.0026
        an=6;
    else if ( (mass >= 12.0) && (mass <= 12.1) )  // Carbon = 12.0110
        an=6;
    else if ( (mass >= 14.0) && (mass <= 14.1) )  // Nitrogen = 14.007
        an=7;
    else if ( (mass >= 15.9) && (mass <= 16.1) )  // Oxygen = 15.9990
        an=8;
    else if ( (mass >= 18.9) && (mass <= 19.0) )  // Fluorine = 18.99800
        an=9;
    else if ( (mass >= 20.0) && (mass <= 20.2) ) // Neon = 20.1797
        an=10;
    else if ( (mass >= 22.9) && (mass <= 23.0) ) // Sodium = 22.989770
        an=11;
    else if ( (mass >= 24.3) && (mass <= 24.4) ) // Magnesium = 24.305000
        an=12;
    else if ( (mass >= 30.0) && (mass <= 31.0) )  // Phosphorus = 30.974000
        an=15;
    else if ( (mass >= 32.0) && (mass <= 32.1) )  // Sulfur = 32.0600
        an=16;
    else if ( (mass >= 35.0) && (mass <= 36.0) )  // Chlorine = 35.45300
        an=17;
    else if ( (mass >= 39.0) && (mass <= 39.2) )  // Potassium = 39.102000
        an=19;
    else if ( (mass >= 40.0) && (mass <= 40.1) )  // Calcium = 40.0800
        an=20;
    else if ( (mass >= 55.0) && (mass <= 56.1) )  // Iron = 55.84700 
        an=26;
    else if ( (mass >= 65.3) && (mass <= 65.4) )  // Zinc = 65.37
        an=30;
    else if ( (mass >= 132.0) && (mass <= 133.0) )  // Cesium = 132.900000
        an=55;
    
    if (an > 0)
        return(an);
    else
        throw(out_of_range("Couldn't identify an atomic number in psf"));
}
