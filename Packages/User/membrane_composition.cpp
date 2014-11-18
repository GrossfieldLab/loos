/*
 *     Compute lipid order params
 *     Assumes the selections are single-bond carbons, and the following
 *     two atoms are hydrogens
 * 
 *     Look at variation of lipid order parameters around a protein
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
#include <boost/format.hpp>
#include <boost/program_options.hpp>


using namespace std;
using namespace loos;
namespace po = boost::program_options;

template<class T> class ValueStore
    {
    private:
        int first, second;
        vector<T> data;
    public:
        ValueStore(const int first_dim, const int second_dim)
            {
            first=first_dim; 
            second=second_dim; 
            if ( (first <= 0) || (second <= 0) )
                {
                throw(runtime_error("dimensions to ValueStore must be >0"));
                }
            data.resize(first*second);
            }

        const T& operator()(const int f, const int s) const
            {
            int index = f*first + s;
            return (data.at(index));
            }

        T& operator()(const int f, const int s)
            {
            int index = f*first + s;
            return (data.at(index));
            }
    };

string system_filename;
string timeseries_filename;
int skip;
string selection;
int first_carbon, last_carbon;
int axis_index;
bool one_res_lipid = false;
bool three_res_lipid = false;

string protein_selection;
int num_x_bins, num_y_bins;
double min_x, max_x, min_y, max_y;

bool dump_timeseries = false;
string traj_filename;

bool block_average = false;
string block_filename;
int ba_first, ba_last;
bool top_leaflet = true;

int frames;

void parseOptions(int argc, char *argv[])
    {
    try
        {
        po::options_description generic("Allowed options");
        generic.add_options()
            ("help,h", "Produce this help message")
            ("timeseries,t", po::value<string>(&timeseries_filename), "File name for outputing timeseries")
            ("x_bins", po::value<int>(&num_x_bins)->default_value(40), "Number of x bins")
            ("y_bins", po::value<int>(&num_y_bins)->default_value(40), "Number of y bins")
            ("min_x", po::value<double>(&min_x)->default_value(-20.0), "Minimum x for histogram")
            ("max_x", po::value<double>(&max_x)->default_value(20.0), "maximum x for histogram")
            ("min_y", po::value<double>(&min_y)->default_value(-20.0), "Minimum x for histogram")
            ("max_y", po::value<double>(&max_y)->default_value(20.0), "maximum x for histogram")
            ("top", po::value<bool>(&top_leaflet)->default_value(true), "top leaflet? (false for bottome leaflet)")
            ;
        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("model", po::value<string>(&system_filename), "Model filename")
            ("traj", po::value<string>(&traj_filename), "Trajectory filename")
            ("sel", po::value<string>(&selection), "Selection string for lipids")
            ("skip", po::value<int>(&skip), "Frames to skip")
            ; 
        po::options_description command_line;
        command_line.add(generic).add(hidden);

        po::positional_options_description p;
        p.add("model", 1);
        p.add("traj", 1);
        p.add("skip", 1);
        p.add("sel", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc,argv).
                options(command_line).positional(p).run(), vm);
        po::notify(vm);

        }
    catch (exception& e)
        {
        cerr << "Error - " << e.what() << endl;
        exit(-1);
        }

    }


int main (int argc, char *argv[])
{

// parse the command line options
parseOptions(argc, argv);

cout << "# " << invocationHeader(argc, argv) << endl;

// Create the data structures for the system and trajectory
AtomicGroup system = createSystem(system_filename);
pTraj traj = createTrajectory(traj_filename, system);

AtomicGroup main_selection = selectAtoms(system, selection);
//split into molecules
vector<AtomicGroup> mol_presplit = main_selection.splitByMolecule();

frames = traj->nframes() - skip;

// skip the equilibration frames
traj->readFrame(skip);

// Allocate space to store the data
double x_bin_width = (max_x - min_x)/num_x_bins;
double y_bin_width = (max_y - min_y)/num_y_bins;

double bin_area = (x_bin_width * y_bin_width);

ValueStore<double> counts(num_x_bins, num_y_bins);

vector<AtomicGroup> molecules;

traj->updateGroupCoords(system);

for (vector<AtomicGroup>::iterator m = mol_presplit.begin(); m!=mol_presplit.end(); ++m)
    {

    if (top_leaflet)
        {
        if ((m->centroid()).z() > 0)
            {
        
            molecules.push_back(*m);

            }
        }
    
    if (!top_leaflet)
        {
        if ((m->centroid()).z() < 0)
            {

            molecules.push_back(*m);

            }
        }
    }



// loop over frames in the trajectory
while (traj->readFrame())
    {
    traj->updateGroupCoords(system);
    
    for (vector<AtomicGroup>::iterator m = molecules.begin(); m!=molecules.end(); ++m)
        {
        double x_val = (m->centroid()).x();
        double y_val = (m->centroid()).y();
        if ( (x_val <= min_x) || (x_val >= max_x) ||
            (y_val <= min_y) || (y_val >= max_y))
            {
            continue;
            }

        int x_bin = (int)((x_val - min_x) / x_bin_width);
        int y_bin = (int)((y_val - min_y) / y_bin_width);
        counts(x_bin, y_bin)++;
        }



    }


cout << "# XBin\tX\tYBin\tY\tCounts" << endl;
for (int k=0; k < num_x_bins; k++)
    {
    float x = min_x + (k+0.5)*x_bin_width;
    for (int m=0; m < num_y_bins; m++)
        {
        float y = min_y + (m+0.5)*y_bin_width;
        cout << k       << "\t"
             << x  << "\t"
             << m       << "\t"
             << y   << "\t"
             << ((counts(k,m) / frames)) / bin_area
             << endl;
        }
    cout << endl;
    }
}

