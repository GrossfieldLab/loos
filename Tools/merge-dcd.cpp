/*
  merge-dcd: combine multiple trajectories into a single long trajectory.  If the target
             trajectory exists, append to it.

  Alan Grossfield
  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School
 
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Alan Grossfield
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

using namespace std;
using namespace loos;

void Usage()
    {
    cerr << "Usage: merge-dcd system recenter-selection output-dcdname "
            "input-dcd [input-dcd2 ...]"
         << endl;    
    cerr << "Giving a empty selection string turns off centering" << endl;
    cerr << "The input dcd files are concatenated in the command line order."
         << endl;
    }

int main(int argc, char *argv[])
    {
    if ( (argc <=1) || 
         ( (argc >=2) && ((strncmp(argv[1], "-h", 2) == 0) ) ) ||
         (argc < 5) )
        {
        Usage();
        exit(-1);
        }

    cout << invocationHeader(argc, argv) << endl;

    AtomicGroup system = createSystem(argv[1]);
    string selection_string = string(argv[2]);
    bool do_recenter = true;
    if ( selection_string.length() == 0 )
        {
        do_recenter = false;
        }

    vector<AtomicGroup> molecules;
    AtomicGroup center;
    vector<AtomicGroup>::iterator m;
    if ( do_recenter )
        {
        center = selectAtoms(system, selection_string);

        if ( system.allHaveProperty(Atom::bondsbit) )
            {
            molecules = system.splitByMolecule();
            }
        else
            {
            molecules = system.splitByUniqueSegid();
            }
        }
    
    DCDWriter output(argv[3], true);
    uint original_num_frames = output.framesWritten();
    cout << "Target trajectory " 
         << argv[3]
         << " has " 
         << original_num_frames
         << " frames."
         << endl;

    uint previous_frames = 0;
    for (int file_index=4; file_index<argc; file_index++ )
        {
        pTraj traj=createTrajectory(argv[file_index], system);
        cout << "File: " << argv[file_index]
             << ": " << traj->nframes();

        if ( previous_frames + traj->nframes() <= original_num_frames) 
            // all of this file is contained in the existing file, skip it
            {
            // increment the frame pointer
            previous_frames += traj->nframes();
            cout << " ( " << previous_frames << " )"
                 << "\tSkipping trajectory " 
                 << endl;

            }
        else
            // we need at least some of the data from this file
            {
            int frames_to_skip = original_num_frames - previous_frames;
            if ( frames_to_skip > 0 )
                {
                traj->seekFrame(frames_to_skip-1);
                }
            else
                {
                frames_to_skip = 0;
                }

            cout << " ( " << previous_frames + traj->nframes() - frames_to_skip
                 << " ) "
                 << "\t Writing " << traj->nframes() - frames_to_skip 
                 << " frames."
                 << endl;

            while ( traj->readFrame() )
                {
                traj->updateGroupCoords(system);

                if ( do_recenter )
                    {
                    GCoord centroid = center.centroid();
                    system.translate(-centroid);
                    for (m=molecules.begin(); m != molecules.end(); ++m )
                        {
                        m->reimage();
                        }
                    }

                output.writeFrame(system);
                previous_frames++;
                }
            }

        }



    }
         
