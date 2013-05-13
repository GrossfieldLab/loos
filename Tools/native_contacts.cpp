/*
  native_contacts: compute the fraction of native contacts in a trajectory,
                   based on an initial structure.
r
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

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage
    {
    public:
        string outfile;
        bool do_output;
        bool exclude_backbone;

        void addGeneric(po::options_description& o) 
            {
            o.add_options()
     ("outfile", po::value<string>(&outfile), "File for timeseries of individual contacts")
     ("exclude-backbone",po::value<bool>(&exclude_backbone)->default_value(false), "Exclude the backbone from contact calculations")
            ;
            }

        bool postConditions(po::variables_map& vm)
            {
            if (vm.count("outfile"))
                {
                do_output=true;
                }
            else
                {
                do_output=false;
                }
            return(true);
            }

    };
// @endcond


string fullHelpMessage(void)
    {
    string s =
"\n"
"    SYNOPSIS\n"
"\n"
"    Report the fraction of native contacts found over the course of \n"
"    a trajectory.\n"
"\n"
"    DESCRIPTION\n"
"\n"
"    The purpose of this tool is to compute the fraction of native contacts\n"
"    found on average over the course of trajectory.  This is intended for\n"
"    use in protein or RNA systems, as a way of tracking the degree to which\n"
"    the molecule is folded.  \n"
"\n"
"    If the model file provided on the command line has coordinates, then \n"
"    those coordinates are used to define \"native\" contacts.  \n"
"    Specifically, the set of atoms to be analyzed is specified on the\n"
"    command line, which is then split by residue.  If the centers of mass\n"
"    of two residues are within the cutoff distance specified on the command\n"
"    line, then those two residues are a native contact.  The same criterion\n"
"    is applied at each successive frame.\n"
"\n"
"    Note: This code does not take periodicity into account.  Rather, it is\n"
"    assumed that the selection (generally a biomolecule such as a protein)\n"
"    is not split across a periodic boundary.  If your simulation code \n"
"    doesn't do this for you, you can use merge-traj with the --fix-imaging\n"
"    flag to set up your trajectory appropriately (although that only works\n"
"    if the selection is a single molecule).\n"
"\n"
"    EXAMPLE\n"
"\n"
"    native_contacts model.psf traj.dcd 5 --selection 'segname == \"PROT\"'\n"
"\n"
"    This uses model.psf as the system file, traj.dcd as the trajectory,\n"
"    sets the cutoff for a native contact at 5 angstroms, and operates on \n"
"    the segment called PROT.  Since PSF files don't have coordinates, the \n"
"    first frame of the trajectory will be used to define which contacts \n"
"    are native.\n"
"\n"
"    If no selection string is provided, then the default is to use "
"    'name == \"CA\"'."
"\n"
"    The output is a time series of the fraction of native contacts present \n"
"    in the trajectory.  At present, there is no option to track the \n"
"    presence of specific contacts (other than using the selection string on\n"
"    the command line to pick out just those contacts).  If you would find\n"
"    this interesting, send email to loos-maintainer@urmc.rochester.edu\n"
"    requesting the feature.\n"
"\n"
        ;
    return(s);
    }

int main (int argc, char *argv[])
{

string header = invocationHeader(argc, argv);

opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
opts::BasicSelection* sopts = new opts::BasicSelection("name == \"CA\"");
opts::RequiredArguments* ropts = new opts::RequiredArguments;
ropts->addArgument("cut", "cutoff");

ToolOptions* topts = new ToolOptions;

opts::AggregateOptions options;
options.add(bopts).add(tropts).add(ropts).add(sopts).add(topts);

if (!options.parse(argc, argv)) 
    {
    exit(-1);
    }

cout << "# " << header << endl;

AtomicGroup system = tropts->model;
pTraj traj = tropts->trajectory;

double cutoff = parseStringAs<double>(ropts->value("cut"));
double cut2 = cutoff*cutoff;


AtomicGroup sel = selectAtoms(system,  sopts->selection);

if (topts->exclude_backbone)
    {
    BackboneSelector backbone;
    NotSelector sidechains(backbone);
    sel = sel.select(sidechains);
    }

vector<AtomicGroup> residues = sel.splitByResidue();

// If they asked for output of individual contacts, set it up
ofstream output;
if (topts->do_output)
    {
    output.open(topts->outfile.c_str());
    if (!output.is_open() ) 
        {
        throw(runtime_error("couldn't open output file"));
        }  
    }

// Check to see if the system has coordinates -- we'll use them if it does,
// otherwise we'll use the first frame of the trajectory as the reference 
// structure.
if ( !(sel[0]->checkProperty(Atom::coordsbit)) )
    {
    traj->readFrame(0);
    traj->updateGroupCoords(system);
    }

// Compute the centers of mass of the selections
uint num_residues = residues.size();
vector<GCoord> centers_of_mass(num_residues);
for (uint i=0; i<num_residues; i++)
    {
    centers_of_mass[i] = residues[i].centerOfMass();
    }

vector<vector<uint> > contacts;
// Find contacts within the threshold distance
for (uint i=0; i<num_residues-1; i++)
    {
    for (uint j=i+1; j< num_residues; j++)
        {
        GCoord diff = centers_of_mass[j] - centers_of_mass[i]; 
        if (diff.length2() <= cut2)
            {
            vector<uint> v(2);
            v[0] = i;
            v[1] = j;
            contacts.push_back(v);
            cout << "# " << (residues[i][0])->resid() << "\t" 
                         << (residues[j][0])->resid() << endl;
            if (topts->do_output)
                {
                output << "# " << (residues[i][0])->resid() << "\t"
                               << (residues[j][0])->resid() << endl;

                }
            }
        }
    }

// Number of native contacts, as a float because we'll need
// to do floating point arithmatic with it later anyway
float num_native_contacts = (float) contacts.size();
cout << "# Total native contacts: " << num_native_contacts << endl;

// Loop over structures in the trajectory
vector<vector<uint> >::iterator p;
int frame = 0;
while (traj->readFrame())
    {
    traj->updateGroupCoords(system);

    // Loop over contacts from the native structure
    int num_contacts = 0;
    for (p=contacts.begin(); p!= contacts.end(); p++)
        {
        uint r1 = p->at(0);
        uint r2 = p->at(1);
        GCoord c1 = residues[r1].centerOfMass();
        GCoord c2 = residues[r2].centerOfMass();
        GCoord diff = c2 - c1;
        if (diff.length2() <= cut2)
            {
            num_contacts++;
            if (topts->do_output) output << "1\t";
            }
        else
            {
            if (topts->do_output) output << "0\t";
            }
        }
    float fraction = num_contacts / num_native_contacts;
    cout << frame << "\t" << fraction << endl;
    if (topts->do_output) output << endl;
    frame++;
    }

}

