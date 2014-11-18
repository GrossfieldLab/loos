/*
 *     Computes bond correlation curves.
 *
 */



/*
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
    cerr << "Usage: bond_correl system traj skip selection "
         << "first_carbon last_carbon max_dT"
         << endl;
    }

int main (int argc, char *argv[])
{
if ( (argc <= 1) || 
     ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
     (argc < 8)
   )
    {
    Usage();
    exit(-1);
    }

cout << "# " << invocationHeader(argc, argv) << endl;

// copy the command line variables to real variable names

char *system_filename = argv[1];
char *traj_filename = argv[2];
int skip = atoi(argv[3]);
char *sel= argv[4];
int first_carbon = atoi(argv[5]);
int last_carbon = atoi(argv[6]);
int max_delta_t = atoi(argv[7]);

// Create the data structures for the system and trajectory
AtomicGroup system = createSystem(system_filename);
pTraj traj = createTrajectory(traj_filename, system);

// NOTE: We assume the selection is a list of all of the relevant carbon atoms. 
//       We'll break it into invidual carbons ourselves (assuming the normal
//       convention of C2, C3, ... 
//       Afterward, we'll figure out the relevant hydrogens ourselves
AtomicGroup main_selection = selectAtoms(system, sel);

// Now break into individual carbons
vector<AtomicGroup> selections;
for (int i =first_carbon; i<=last_carbon; i++)
    {
    string sel_string = string(sel);
    char carbon_name[4];
    sprintf(carbon_name, "%d", i);
    string name = string(" && name == \"C") + string(carbon_name)
                  + string("\"");
    sel_string.insert(sel_string.size(), name);
    //Parser p(sel_string.c_str());
    //KernelSelector parsed(p.kernel());
    //selections.push_back(main_selection.select(parsed));
    selections.push_back(selectAtoms(main_selection, sel_string.c_str()));
    }

// Now, figure out which hydrogens go with each carbon selected
// hydrogen_list mimics the structure of selections, so the hydrogens
// bound to the jth carbon of the ith selection will be found at
// hydrogen_list[i][j].
vector< vector<AtomicGroup> > hydrogen_list(selections.size());
HydrogenSelector hyd_sel;


// bond_counts stores the number of bonds for a given carbon,
// so we know how many to increment through when filling the
// larger vector
vector<int> bond_counts;

for (unsigned int i=0; i<selections.size(); i++) 
    {
    AtomicGroup *s = &(selections[i]);
    AtomicGroup::Iterator iter(*s);
    pAtom p;
    int count = 0;
    while (p = iter()) 
        {
        vector<int> atom_ids = p->getBonds();
        AtomicGroup bonded = system.groupFromID(atom_ids);
        AtomicGroup bonded_hydrogens = bonded.select(hyd_sel);
        hydrogen_list[i].push_back(bonded_hydrogens);
        count = count + bonded_hydrogens.size();
        }
    bond_counts.push_back(count);
    }


// skip the equilibration frames
traj->readFrame(skip);

// Initialize the data vector that will store timeseries
// information for each carbon-hydrogen bond
vector< vector< vector<GCoord> > > series;

for (unsigned int i = 0; i < selections.size(); i ++)
    {

    vector < vector < GCoord > > temp_carbon;

    for (int j = 0; j < bond_counts.at(i); j++)
        {

        vector<GCoord> temp_bond;
        temp_carbon.push_back(temp_bond);

        }

    series.push_back(temp_carbon);
    
    }

// loop over pdb files
while (traj->readFrame())
    {
    traj->updateGroupCoords(system);

    // loop over sets of selected carbons
    for (unsigned int i=0; i<selections.size(); i++)
        {

        AtomicGroup *g = &(selections[i]);
        
        // keeps track of bond within this selection to make sure
        // we store the data appropriately
        int bond_instance = 0;
        for (int j=0; j < g->size(); j++)
            {
            // get the carbon
            pAtom carbon = g->getAtom(j);
            // get the relevant hydrogens
            AtomicGroup *hyds = &(hydrogen_list[i][j]);
            
            AtomicGroup::Iterator iter(*hyds);
            pAtom h;
            while (h = iter() )
                {
                GCoord v = carbon->coords() - h->coords();
                GCoord prod = v / v.length();
                ( ( series.at(i) ).at(bond_instance) ).push_back(prod);
                bond_instance++;
                }
            }
        }
        
    }

// Print out carbon names
cout << "#dT";

// First line of output with column titles (carbon names)
for (unsigned int i = 0; i < selections.size(); i++)
    {
    pAtom pa = selections[i].getAtom(0);
    cout << "\t" << pa->name() << "\t" << pa->name() << "dev";
    }

    cout << endl;

// Loop over all t increments
for (int l = 0; l <= max_delta_t; l++)
    {
    cout << l;
    // Loop over all carbons
    for (unsigned int i = 0; i < series.size(); i++)
        {

        // Vector of the averages for each bond
        vector<double> bond_set;

        // Loop over all bonds for each carbon
        for (unsigned int j = 0; j < (series.at(i)).size(); j++)
            {
            // Store a counter and average for each bond
            double average_Y20 = 0.0;

            double average_real_Y21 = 0.0;
            double average_imag_Y21 = 0.0;
            double average_Y21 = 0.0;

            double average_real_Y22 = 0.0;
            double average_imag_Y22 = 0.0;
            double average_Y22 = 0.0;
            
            int counter = 0;

            // calculate average vector values for subtraction
            for (unsigned int k = 0; k < ((series.at(i)).at(j)).size(); k++)
                {
                double z = (((series.at(i)).at(j)).at(k)).z();
                double y = (((series.at(i)).at(j)).at(k)).y();
                double x = (((series.at(i)).at(j)).at(k)).x();
                
                // Calculate Y20 average
                average_Y20 =+ (((3/2) * z * z) - 0.5) * (((3/2) * z * z) - 0.5);
                
                // Caculate the average for Y21. This involves calculating the average
                // of the real and imaginary components separately. Then we take the magnitude
                // of this (sqrt(x*x + y*y)). Then we square this magnitude. The magnitude and
                // the square of this can be combined to just be (x*x + y*y).
                //
                // The extra coefficient is ignored until later, since it drops out anyway.
                
                // Sum real and imaginary components
                average_real_Y21 =+ z * sqrt(1-(z*z)) * x;
                average_imag_Y21 =+ z * sqrt(1-(z*z)) * y;

                // Y22 is calculated the same way as Y21, though the actual equation is different.
                average_real_Y22 =+ (1-(z*z))*((2*z*z)-1);
                average_imag_Y22 =+ (1-(z*z))*(2*z*(sqrt(1-(z*z))));

                counter++;
                }

            average_Y20 = average_Y20 / counter;
            
            // Average of real and imag components for Y21
            average_real_Y21 = average_real_Y21 / counter;
            average_imag_Y21 = average_imag_Y21 / counter;

            // Average of real and imag components for Y22
            average_real_Y22 = average_real_Y22 / counter;
            average_imag_Y22 = average_imag_Y22 / counter;

            // Square of magnitude for Y21
            // Includes 1.5 factor, which is simply the coefficients squared, multiplied by the 4pi/5 in the
            // main equation (drops out to a simple value).
            // Just includes the real component for now
            average_Y21 = 1.5* ((average_real_Y21 * average_real_Y21) - (average_imag_Y21 * average_imag_Y21));

            // Square of magnitude for Y22
            // Just includes real component for now
            average_Y22 = (3/8)*((average_real_Y22 * average_real_Y22) - (average_imag_Y22 * average_imag_Y22));

            // Store a counter and average for each bond
            double bond_final = 0.0;
            int bond_counter = 0;
           
            double Y20 = 0.0;
            double Y21 = 0.0;
            double Y22 = 0.0;


            double average_dt_Y20 = 0.0;
            double average_dt_Y21 = 0.0;
            double average_dt_Y22 = 0.0;

            // Loop over the timeseries for each bond
            for (int k = 0; k < ((int)((series.at(i)).at(j)).size()) - l; k++)
                {
                double z_t = (((series.at(i)).at(j)).at(k)).z();
                double z_dt = (((series.at(i)).at(j)).at(k+l)).z();
                double y_t = (((series.at(i)).at(j)).at(k)).y();
                double y_dt = (((series.at(i)).at(j)).at(k+l)).y();
                double x_t = (((series.at(i)).at(j)).at(k)).z();
                double x_dt = (((series.at(i)).at(j)).at(k)).z();

                double y_t2 = y_t*y_t;
                double y_dt2 = y_dt*y_dt;

                double z_t2 = z_t*z_t;
                double z_dt2 = z_dt*z_dt;

                // Calculate Y20
                Y20 = ((1.5*z_t2) - 0.5)*((1.5*z_dt2) - 0.5);
                average_dt_Y20 += Y20;

                // Calculate Y21
                double z_components = 1.5 * (sqrt(1-z_t2) * z_t) * (sqrt(1-z_dt2) * z_dt);
                double real_dt_Y21 = z_components * ((x_t*x_dt) - (y_t*y_dt));
                double imag_dt_Y21 = z_components * ((y_t*x_dt) - (x_t*y_dt));
                //Y21 = real_dt_Y21;
                Y21 = (real_dt_Y21 * real_dt_Y21) - (imag_dt_Y21 * imag_dt_Y21);

                // Calculate Y22
                z_components = (3/8) * (1-z_t2) * (1-z_dt2);
                double real_dt_Y22 = z_components * (((1-(2*y_t2))*(1-(2*y_dt2))) + ((2*y_t*x_t)*(2*y_dt*x_dt)));
                double imag_dt_Y22 = z_components * (((2*y_t*x_t)*(1-(2*y_dt2))) - ((1-(2*y_t2))*(2*y_dt*x_dt)));
                Y22 = (real_dt_Y22 * real_dt_Y22) - (imag_dt_Y22 * imag_dt_Y22);
                //Y22 = real_dt_Y22;
                bond_counter++;    

                }

            average_dt_Y20 = average_dt_Y20 / bond_counter;
            double final_Y20 = average_dt_Y20 - average_Y20;

            average_dt_Y21 = average_dt_Y21 / bond_counter;
            double final_Y21 = average_dt_Y21 - average_Y21;

            average_dt_Y22 = average_dt_Y22 / bond_counter;
            double final_Y22 = average_dt_Y22 - average_Y22;

            bond_final = (final_Y20 + (2*final_Y21) + (2*final_Y22)) / 5;

            bond_set.push_back(bond_final);

            }
        // Calculate average and output stdev
        TimeSeries<double> t_set;
        t_set = TimeSeries<double>(bond_set);

        cout << "\t" << t_set.average() << "\t" 
            << t_set.stdev();
        }
    cout << endl;
    }

}
