/* 
  Locate waters which cross the membrane

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


#include <iostream>
#include <map>
#include <loos.hpp>

using namespace std;


void Usage()
    {
    cerr << "Usage: crossing_waters system traj inner_threshold outer_threshold"
         << endl;
    }

class InternalWater
    {
    public:
        InternalWater(pAtom at, int first, double z) 
            : atom(at),
              entry_frame(first),
              exit_frame(-1),
              exited_to_positive(false)
              {
              set_entered_side(z);
              }
    
  
        int entered() const {return(entry_frame);}

        void set_entered_side(double z)
            {
            if (z>0)
                {
                entered_from_positive = true;
                }
            else
                {
                entered_from_positive = false;
                }
            }

        void set_exit_side(double z)
            {
            if (z>0)
                {
                exited_to_positive = true;
                }
            else
                {
                exited_to_positive= false;
                }
            }

        bool from_pos() const {return(entered_from_positive);}
        int exited() const {return(exit_frame);}
        void exit(const int frame) {exit_frame = frame;}

        bool crossed()
            {
            if (exit_frame > 0)  // we have exited the membrane
                {
                if ( (entered_from_positive  && !exited_to_positive) ||
                     (!entered_from_positive &&  exited_to_positive) )
                    {
                    return(true);
                    }
                }
            return(false);
            }

        int lifetime()
            {
            if (exit_frame > 0)
                {
                return( exited() - entered() );
                }
            return(-1);
            }

        pAtom get_atom()
            {
            return(atom);
            }

        int atom_id()
            {
            return(atom->id());
            }

    private:
        pAtom atom;
        int entry_frame;
        int exit_frame;
        bool entered_from_positive; 
        bool exited_to_positive;
    };

int main (int argc, char *argv[])
{
if ( (argc <= 1) ||
     ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
     (argc < 5)
   )
    {
    Usage();
    exit(-1);
    }

cout << "# " << loos::invocationHeader(argc, argv) << endl;

AtomicGroup system = loos::createSystem(argv[1]);
pTraj traj = loos::createTrajectory(argv[2], system);
greal inner_threshold = atof(argv[3]);
greal outer_threshold = atof(argv[4]);

// Select the water oxygens
HeavySolventSelector water_heavy; // select atoms which are both water 
                                  // and not hydrogens
AtomicGroup water = system.select(water_heavy); // apply the selection


// Set up a map to store the waters
map<pAtom,InternalWater> internal_waters;
vector<InternalWater> exited_waters;

int frame = 0;
while (traj->readFrame())
    {
    traj->updateGroupCoords(system);
    bool inside_inner, inside_outer;

    // loop over water oxygens
    pAtom w;
    AtomicGroup::Iterator iter(water);
    while (w=iter())
        {
        double z = w->coords().z();
        if (fabs(z) < inner_threshold)
            {
            inside_inner = true;
            inside_outer = true;
            }
        else if (fabs(z) < outer_threshold)
            {
            inside_inner = false;
            inside_outer = true;
            }
        else
            {
            inside_inner = false;
            inside_outer = false;
            }
        
        // is this a water which is currently inside?
        map<pAtom,InternalWater>::iterator wat = internal_waters.find(w);
        if (wat != internal_waters.end())
            {
            // Yes, the water is listed as internal
            if (inside_outer)
                {
                // do nothing -- it's still inside
                }
            else
                {
                // the water has left -- mark it so and clean up
                InternalWater &i = wat->second;
                i.exit(frame);
                i.set_exit_side(z);
                exited_waters.push_back(i);
                internal_waters.erase(w);
                // TODO: debugging code, don't leave in
                //if (internal_waters.find(w) != internal_waters.end())
                //    {
                //    throw(runtime_error("Erase didn't erase"));
                //    }
                }
            }
        else
            // No, the water wasn't previously inside
            {
            if (inside_inner)
                {
                // create a new InternalWater and store it
                internal_waters.insert(
                    pair<pAtom,InternalWater>(w, InternalWater(w, frame, z)));
                }
            else
                {
                // not inside: do nothing
                }
            }
        }
        frame++;
    }

vector<InternalWater>::iterator wat;
cout << "# Total frames = " << frame << endl;
cout << "#AtomID\tLifetime\tEntered\tExited" << endl;
for (wat = exited_waters.begin(); wat != exited_waters.end(); wat++)
    {
    if (wat->crossed())
        {
        cout << wat->atom_id() << "\t"
             << wat->lifetime() << "\t"
             << wat->entered() << "\t"
             << wat->exited() 
             << endl;
        }
    }

}
