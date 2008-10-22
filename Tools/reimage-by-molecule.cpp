/*
   reimage-by-molecule
  
   Read a psf and dcd, reimage each frame by molecule, and write a new
   dcd

   Alan Grossfield
 
  
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
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

#include <loos.hpp>

void Usage()
    {
    cerr << "Usage: reimage-by-molecule psf indcd outdcd [xbox ybox zbox]"
         << endl;
    } 

int main(int argc, char *argv[]) 
    {
    if ( (argc <= 1) ||
         ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
         (argc < 4)
       )
        {
        Usage();
        exit(-1);
        }

    cout << "# " << invocationHeader(argc, argv) << endl;

    PSF psf(argv[1]);
    DCD dcd(argv[2]);
    DCDWriter dcd_out(argv[3]);

    if (!dcd.hasPeriodicBox())
        {
        // if the dcd doesn't have periodic box info, we need it from
        // the user
        if (argc < 7)
            {
            cerr << "DCD " << argv[2] << " doesn't have box info" << endl;
            cerr << "so you'll need to supply it on the command line" << endl;
            exit(-1);
            }

        // use the supplied box, set it in psf
        double xbox = atof(argv[4]);
        double ybox = atof(argv[5]);
        double zbox = atof(argv[6]);
        GCoord box(xbox, ybox, zbox);
        psf.periodicBox(box);
        }

    // duplicate the info from dcd in dcd_out
    dcd_out.setHeader(dcd.natoms(), dcd.nsteps(), dcd.timestep(), 
                     true );
    dcd_out.setTitles(dcd.titles());
    dcd_out.writeHeader();

    // split the system by molecule
    vector<AtomicGroup> molecules = psf.splitByMolecule();

    // split the system by segid
    vector<AtomicGroup> segments = psf.splitByUniqueSegid();


    // Loop over the frames of the dcd and reimage each molecule
    vector<AtomicGroup>::iterator m;
    while (dcd.readFrame())
        {
        dcd.updateGroupCoords(psf);
        for (m = segments.begin(); m != segments.end(); m++)
            {
            m->reimage();
            }
        for (m = molecules.begin(); m != molecules.end(); m++)
            {
            m->reimage();
            }
        dcd_out.writeFrame(psf);
        }

    }

