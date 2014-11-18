/*
 *Map the average phosphate Z value (assumes protein or bilayer is centered 
 *at zero to be useful). If groups larger than a single atom are given, the
 *returned data will be for the centroid of the molecules or residues.
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
int axis_index;

string protein_selection;
int num_x_bins, num_y_bins;
double min_x, max_x, min_y, max_y;
string split;


string traj_filename;

void parseOptions(int argc, char *argv[])
    {
    try
        {
        po::options_description generic("Allowed options");
        generic.add_options()
            ("help,h", "Produce this help message")
            ("prot_select,p", po::value<string>(&protein_selection)->default_value(string("name == \"CA\")")), "Selection of atoms defining the protein")
            ("x_bins", po::value<int>(&num_x_bins)->default_value(40), "Number of x bins")
            ("y_bins", po::value<int>(&num_y_bins)->default_value(40), "Number of y bins")
            ("min_x", po::value<double>(&min_x)->default_value(-40.0), "Minimum x for histogram")
            ("max_x", po::value<double>(&max_x)->default_value(40.0), "maximum x for histogram")
            ("min_y", po::value<double>(&min_y)->default_value(-40.0), "Minimum x for histogram")
            ("max_y", po::value<double>(&max_y)->default_value(40.0), "maximum x for histogram")
            ("split_by", po::value<string>(&split)->default_value("by-molecule"), "how to split the targets (by-molecule or by-residue)");
            ;
        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("model", po::value<string>(&system_filename), "Model filename")
            ("traj", po::value<string>(&traj_filename), "Trajectory filename")
            ("sel", po::value<string>(&selection), "Selection for which to calculate density")
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

        if (vm.count("help")              ||
                !vm.count("model")        || !vm.count("traj")         ||
                !vm.count("skip")         || !vm.count("sel")          
           )
            {
            cerr << "Usage: " << argv[0] << " "
                 << "model-name trajectory-name skip-frames selection-string "
                 << endl;
            cerr << generic;
            exit(-1);
            }

        }
    catch(exception& e)
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

AtomicGroup protein = selectAtoms(system, protein_selection);
//AtomicGroup reference = protein.copy();

AtomicGroup target = selectAtoms(system, selection);

vector<AtomicGroup> targets = target.splitByMolecule();

if (split=="by-residue"){

    targets = target.splitByResidue();

}

//GCoord center = reference.centroid();
//reference.translate(-center);



// skip the equilibration frames
traj->readFrame(skip);

const int num_frames = traj->nframes() - skip;

// Allocate space to store the data
double x_bin_width = (max_x - min_x)/num_x_bins;
double y_bin_width = (max_y - min_y)/num_y_bins;

ValueStore<int> counts_upper(num_x_bins, num_y_bins);
ValueStore<int> counts_lower(num_x_bins, num_y_bins);
ValueStore<double> sum_upper(num_x_bins, num_y_bins);
ValueStore<double> sum_lower(num_x_bins, num_y_bins);
// For normalization purposes, I'll track the total counts
// for each frame
int total_counts = 0;

// loop over frames in the trajectory
int frame_index = 0;
while (traj->readFrame())
    {
    traj->updateGroupCoords(system);

    // rotate the system so that the protein is aligned against reference
    // TODO: not perfect, because this could tilt the protein -- not sure how to
    //       handle it.
//    GMatrix matrix = protein.superposition(reference);
//    XForm transform(matrix);
//    system.applyTransform(transform);
//    system.reimageByAtom();
    GCoord prot_centroid = protein.centroid();
    prot_centroid.z() = 0.0;

    // loop over all molecules in selection carbons
    for (int i=0; i<targets.size(); i++)
        {
        
        GCoord target_vector = targets[i].centroid() - prot_centroid;
        if ( (target_vector.x() <= min_x) || (target_vector.x() >= max_x) ||
             (target_vector.y() <= min_y) || (target_vector.y() >= max_y) )
            {
            continue;
            }

        int x_bin = (int)((target_vector.x() - min_x) / x_bin_width);
        int y_bin = (int)((target_vector.y() - min_y) / y_bin_width);

        double diff = targets[i].centroid().z() - target.centroid().z();
        if ( diff > 0)
            {
            counts_upper(x_bin, y_bin)++;
            sum_upper(x_bin, y_bin) = sum_upper(x_bin, y_bin) + diff;
            total_counts++;   
            } else {
            counts_lower(x_bin, y_bin)++;
            sum_lower(x_bin, y_bin) = sum_lower(x_bin, y_bin) + diff;
            total_counts++;
            }
        }
    }

float bin_area = (x_bin_width * y_bin_width);

cout << "# XBin\tX\tYBin\tY\tUpper\tLower" << endl;

double average = ((double)total_counts / ((max_x - min_x) * (max_y - min_y))) / (double)num_frames;

for (int k=0; k < num_x_bins; k++)
    {
    float x = min_x + (k+0.5)*x_bin_width;
    for (int m=0; m < num_y_bins; m++)
        {
        float y = min_y + (m+0.5)*y_bin_width;

        double sum_upper_normalized = 0;
        double sum_lower_normalized = 0;
        
        if (counts_upper(k,m) > 0)
            {
            sum_upper_normalized = (((double)sum_upper(k,m)) / (double)counts_upper(k,m));
            }

        if (counts_lower(k,m) > 0)
            {
            sum_lower_normalized = (((double)sum_lower(k,m)) / (double)counts_lower(k,m));
            }

        double sum_total = sum_upper_normalized + abs(sum_lower_normalized);

        cout << k       << "\t"
             << x  << "\t"
             << m       << "\t"
             << y   << "\t"
             << sum_upper_normalized << "\t"
             << sum_lower_normalized << "\t"
             << sum_total << "\t"
             << endl;
        }
    cout << endl;
    }
}

