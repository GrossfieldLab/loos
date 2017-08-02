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
using namespace loos;



void Usage()
    {
    cerr << "Usage: crossing_waters system traj inner_threshold outer_threshold [water-selection]"
         << endl;
    }

// @cond TOOLS_INTERNAL

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

        const int exit_side() const
            {
            if (exited_to_positive)
                {
                return +1;
                }
            else
                {
                return -1;
                }
            }

    private:
        pAtom atom;
        int entry_frame;
        int exit_frame;
        bool entered_from_positive;
        bool exited_to_positive;
    };

// @endcond

string fullHelpMessage(void)
    {
    string s =
"\n"
"SYNOPSIS\n"
"\n"
"Track the rate at which water molecules cross the membrane.\n"
"\n"
"DESCRIPTION\n"
"\n"
"This tool measures the rate at which water molecules cross the lipid membrane.  \n"
"To simplify matters, the tool assumes that the membrane is centered at z=0, \n"
"and keeps track of water molecules passing through z=0.  To differentiate \n"
"between waters passing through the membrane and those that simply pass through \n"
"the periodic boundary, the user specifies 2 values, inner_threshold and \n"
"outer_threshold, which specify a distance from the membrane center at which \n"
"the waters are determined to have entered and exited the membrane, \n"
"respectively.  The rationale for using 2 thresholds is that we only want to \n"
"track waters with a reasonable chance of crossing the membrane, so we use a \n"
"restrictive threshold there (the water has to really be in the membrane before \n"
"we pay attention), but on the other hand we don't want to say it's out until \n"
"it's safely outside the membrane, so we use a larger threshold there.  The \n"
"optimum choice for these values depends on the thickness of the membrane, but \n"
"10 and 20 are a reasonable start.\n"
"\n"
"There's an additional optional flag to control how water is selected, in\n"
"case your naming conventions are different from ours.  Just add a selection\n"
"string after the outer_threshold.  The code will internally remove any\n"
"hydrogens, so you don't have to put that in your selection unless you want.\n"
"\n"
"The output is a list of waters that crossed the membrane, how many frames \n"
"each water spent inside the membrane, and the frames it entered and left.  The \n"
"final column is either 1 or -1; the former indicated that the water exited on \n"
"the +z side of the membrane, the latter the -z side.\n"
"\n"
"EXAMPLE\n"
"\n"
"crossing-waters system.psf traj.dcd 10.0 20.0\n"
"\n"
"This will read system.psf and the trajectory file traj.dcd, and use 10 and 20 \n"
"angstroms as the inner and outer threshold.  The output will look like this:\n"
"\n"
"# crossing-waters 'system.psf' 'traj.dcd' '10' '20' - alan (Tue Mar 13 14:32:49 2012) {/directory/you/were/working/in} [2.0.0 120313]\n"
"# Total frames = 719\n"
"#AtomID	Lifetime	Entered	Exited	ExitedPositive\n"
"38546	2	0	2	-1\n"
"25136	2	4	6	1\n"
"25856	1	5	6	-1\n"
"35909	1	8	9	-1\n"
"39665	1	8	9	-1\n"
"\n"
        ;
    return(s);
    }


int main (int argc, char *argv[])
{
if ( (argc >=2) )
    {
    if (string(argv[1]) == string("-h"))
        {
        Usage();
        exit(-1);
        }
    else if (string(argv[1]) == string("--fullhelp"))
        {
        cerr << fullHelpMessage() << endl;
        exit(-1);
        }
    }
if (argc < 5)
    {
    Usage();
    exit(-1);
    }

cout << "# " << invocationHeader(argc, argv) << endl;

AtomicGroup system = createSystem(argv[1]);
pTraj traj = createTrajectory(argv[2], system);
greal inner_threshold = atof(argv[3]);
greal outer_threshold = atof(argv[4]);

bool manual_water = false;
char *water_selection;
if (argc > 5)
    {
    manual_water = true;
    water_selection = argv[5];
    }

// Select the water oxygens
HeavySolventSelector water_heavy; // select atoms which are both water
                                  // and not hydrogens

AtomicGroup water;
if (manual_water)  // apply the supplied selection, but ensure no hydrogens
    {
    water = selectAtoms(system, water_selection);
    water = selectAtoms(water, "!hydrogen");
    }
else
    {
    water = system.select(water_heavy); // apply the selection
    }


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
    while ( (w=iter()) )
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
cout << "# Number of waters = " << water.size() << endl;
cout << "#AtomID\tLifetime\tEntered\tExited\tExitedPositive" << endl;
for (wat = exited_waters.begin(); wat != exited_waters.end(); wat++)
    {
    if (wat->crossed())
        {
        cout << wat->atom_id() << "\t"
             << wat->lifetime() << "\t"
             << wat->entered() << "\t"
             << wat->exited() << "\t"
             << wat->exit_side()
             << endl;
        }
    }

}
