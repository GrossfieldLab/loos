/*
  micelle
  
  (c) 2009 Joshua N. Horn, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry

  Determines molecule contacts to calculate the number of micelles
  formed and other properties of these micelles.
*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008-2009, Tod D. Romo
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

bool emptyAtomicGroup (AtomicGroup g) { return (g.size() == 0); }

bool contact(const AtomicGroup& A, const AtomicGroup& B, const GCoord& box)
    {
;
    bool value = false;

    for (AtomicGroup::const_iterator j = A.begin(); (j != A.end()) && (value != true); ++j)
        {

        GCoord u = (*j)->coords();

        for (AtomicGroup::const_iterator k = B.begin(); (k != B.end()) && (value != true); ++k)
            {

            GCoord v = (*k)->coords();

            if (u.distance(v, box) < 7)
                {

                value = true;

                }
                

            }

        }

    return(value);

    }

int main (int argc, char *argv[])
{
// Build usage if/then here



// Print the command line arguments
cout << "# " << invocationHeader(argc, argv) << endl;

// Create the system and the trajectory file
AtomicGroup system = createSystem(argv[1]);
pTraj traj = createTrajectory(argv[2], system);

// String describing the first selection
char *selection = argv[3];

// Get the number of frames to discard as equilibration
int skip = atoi(argv[4]);

// Get the last frame to analyze
int lastFrame = atoi(argv[5]);
if (lastFrame == 0)
    {

        lastFrame = traj->nframes();

    }

vector<AtomicGroup> molecules;

molecules = system.splitByResidue();

//cout << molecules.size() << endl;

cout << "# frame aggregates (mols in each)" << endl;

// Set up the selector to define the selected group
Parser parser(selection);
KernelSelector parsed_sel(parser.kernel());

// Loop over the molecules and add them to selection
vector<AtomicGroup> molecule_groups;
vector<AtomicGroup>::iterator m;

// Build atomic groups and push them into molecule_groups
for (m=molecules.begin(); m!=molecules.end(); m++)
    {
    AtomicGroup tmp = m->select(parsed_sel);
    if (tmp.size() > 0)
        {
        molecule_groups.push_back(tmp);

        }

    }

//cout << "size: " << molecule_groups.size() << endl;

// Skip the initial frames as equilibration
traj->readFrame(skip); 

// read the initial coordinates into the system
traj->updateGroupCoords(system);

// "Connectivity tree" for micelles

int framecount = skip;
while (traj->readFrame() && framecount < lastFrame)
    {
    cout << framecount;

    // "Connectivity tree" for micelles
    vector< vector<int> > mols;
    for (unsigned int p = 0; p<molecule_groups.size(); p++)
        {

        vector<int> tmp;
        mols.push_back(tmp);

        }

    traj->updateGroupCoords(system);

    // loop over groups and find contacts
    vector<AtomicGroup>::iterator l;
    for (unsigned int l=0; l!=molecule_groups.size(); l++)
        {

        vector<AtomicGroup>::iterator n;
        for (unsigned int n=l + 1; n!=molecule_groups.size(); n++)
            {

            bool in_contact = false;

            // to speed up program, do a quick test on center of mass, only
            // calculate contact if within some threshold
            //
            GCoord tmp_box = system.periodicBox();

                if (contact(molecule_groups.at(n), molecule_groups.at(l), tmp_box) == true)
                    {
                    in_contact = true;
                    }


            if (in_contact == true)
                {

                (mols.at(l)).push_back(n);
                

                }

            }

        }


        //for (int p = 0; p < mols.size(); p++)
          //  {

            //cout << p << " [";
            //for (int q = 0; q < (mols.at(p)).size(); q++)
                //{

                //cout << (mols.at(p)).at(q) << " ";

              //  }
            //cout << "] " << endl;

            //}

        // build micelles
        vector<AtomicGroup> molecules;
        vector<AtomicGroup> micelles;

        // fill micelles with the micelles of 
        for (unsigned int a = 0; a < mols.size(); a++)
            {

            AtomicGroup temp = molecule_groups[a];
            for (unsigned int b = 0; b < (mols.at(a)).size(); b++)
                {

                temp = temp.merge(molecule_groups[(mols.at(a)).at(b)]);

                }

            molecules.push_back(temp);

            }

        while (molecules.empty() == false)
            {
            //cout << micelles.size() << " ";
            //cout << molecules.size() << " " << endl;

            bool changed = true;
            while (changed == true)
                {

                changed = false;
                for (unsigned int i = 1; i < molecules.size(); i++)
                    {

                    if (((molecules.at(0)).intersect(molecules.at(i))).size() != 0)
                        {
                        changed = true;
                        //cout << "combining 0 " << i << endl;
                        molecules.at(0) = ((molecules.at(0)).merge(molecules.at(i)));
                        AtomicGroup emptyGroup;
                        molecules.at(i) = emptyGroup;
                        }

                    }

                }
                micelles.push_back(molecules.at(0));
                AtomicGroup emptyGroup;
                molecules.at(0) = emptyGroup;
           
                //cout << "cleanup" << endl;
            // cleanup!
                
                vector<AtomicGroup>::iterator toDelete = remove_if(molecules.begin(), molecules.end(), emptyAtomicGroup);
                molecules.erase(toDelete, molecules.end());
            }            

        int total = 0;
        cout << " " << micelles.size();
        for(unsigned int t = 0; t < micelles.size(); t++)
            {

            total = total + ((micelles.at(t)).size() / 4);
            cout << " " << ((micelles.at(t)).size() / 4);
            //cout << " [" << ((micelles.at(t)).getAtom(0))->id() << "] ";

            }
        cout << endl;
       // cout << " total: " << total << endl;
    
        framecount++;
    } // end frame while

} // end main

