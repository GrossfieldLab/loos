/*
 *   psf.cpp
 *   (c) 2008 Alan Grossfield
 *   Department of Biochemistry and Biophysics
 *   University of Rochester Medical School
 *
 *   Simple CHARMM/NAMD PSF file reader
 *
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

    // Is the atom fixed or mobile?
    // for now, we're going to silently drop this
    ss >> fixed;  

    // now intialize the rest of the stuff
    pa->recordName(string(""));
    GCoord c= GCoord(99999.99, 99999.99, 99999.99); // Marks
						    // uninitialized
						    // coords 
    pa->coords(c);

    append(pa);
}
