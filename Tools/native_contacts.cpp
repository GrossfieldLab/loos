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

#include "loos.hpp"

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage
    {
    public:
        string outfile;
        string reference;
        string per_residue_filename;
        bool do_output;
        bool exclude_backbone;
        bool use_periodicity;
        bool use_reference;
        bool do_per_residue;
        bool exclude_consecutive;

        void addGeneric(po::options_description& o)
            {
            o.add_options()
     ("outfile", po::value<string>(&outfile), "File for timeseries of individual contacts")
     ("exclude-backbone",po::value<bool>(&exclude_backbone)->default_value(false), "Exclude the backbone from contact calculations")
     ("periodic", "Use periodicity when computing contacts")
     ("reference", po::value<string>(&reference), "Coordinate file to use as reference structure")
     ("per-residue", po::value<string>(&per_residue_filename), "Output per-residue native contact frequency to this file")
     ("exclude-consecutive", "Exclude consecutive residues")
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

            if (vm.count("periodic"))
                {
                use_periodicity = true;
                }
            else
                {
                use_periodicity = false;
                }

            if (vm.count("exclude-consecutive"))
                {
                exclude_consecutive = true;
                }
            else
                {
                exclude_consecutive = false;
                }

            if (vm.count("reference"))
                {
                use_reference = true;
                }
            else
                {
                use_reference = false;
                }

            if (vm.count("per-residue"))
                {
                do_per_residue = true;
                }
            else
                {
                do_per_residue = false;
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
"    By default, the model file provided on the command line has coordinates, then \n"
"    those coordinates are used to define \"native\" contacts.  \n"
"    If the model file doesn't have coordinates, then the first frame of the\n"
"    trajectory is used.\n"
"\n"
"    Alternatively, you can supply a separate structure containing reference\n"
"    coordinates (e.g. a pdb file with the original crystal coordinates).\n"
"    The only restriction is that the same selection string that picks out\n"
"    the residues of interest from the system file must also apply to the \n"
"    reference file.\n"
"\n"
"    The set of atoms to be analyzed is specified on the\n"
"    command line, which is then split by residue.  If the centers of mass\n"
"    of two residues are within the cutoff distance specified on the command\n"
"    line, then those two residues are a native contact.  The same criterion\n"
"    is applied at each successive frame.\n"
"\n"
"    Note: This code does not take periodicity into account by default,\n"
"    because in most cases (e.g. a protein or RNA) the molecule will be \n"
"    in a single unit cell.  If you want periodicity, add the flag \n"
"    '--periodic' on the command line.  If you give this flag and supply an \n"
"    initial structure that does not have box information, you will get a \n"
"    warning, and the initial identification of contacts will be done without\n"
"    using the periodic image.  If this is not the desired behavior, you'll \n"
"    need to add the box information to the initial structure by hand first,\n"
"    or use the first frame of the trajectory as the reference.\n"
"\n"
"    The --exclude-consecutive option causes the code to ignore residues\n"
"    consecutive in sequence when computing the list of native contacts.\n"
"    Note: this is done in a naive way, without checking that the consecutive\n"
"    residues are part of the same chain.  "
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
"    If no selection string is provided, then the default is to use\n "
"    'name == \"CA\"'."
"\n"
"    In addition, one can select just the sidechains using the\n "
"    --exclude-backbone flag; this can be combined with other selections.\n"
"    Turn it one with --exclude-backbone 1 \n"
"\n"
"    If you supply the \"--outfile\" option, you will also get a time series for \n"
"    all of the individual pairs of residues.\n"
"\n"
"    If you supply the \"--per-residue FILENAME\", the program will output \n"
"    the average fractional native contacts for each residue to FILENAME.\n"
"    Residues with no native contacts will have a value of -1.\n"
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

// Figure out what to use as a reference structure
// --If one was supplied on the command line use that, and
//   use the same selections applied to the main system to
//   select out the equivalent group.
// --If there was no reference structure, check to see if the
//   model file has coordinates, and use those.
// --If neither is true, use frame 0 of the trajectory.
if (topts->use_reference)
    {
    AtomicGroup reference = createSystem(topts->reference);
    AtomicGroup ref_sel = selectAtoms(reference, sopts->selection);
    if (topts->exclude_backbone)
        {
        BackboneSelector backbone;
        NotSelector sidechains(backbone);
        ref_sel = ref_sel.select(sidechains);
        }
    // copy the coordinates from ref_sel to sel, after sanity checking
    if (ref_sel.size() != sel.size())
        {
        cerr << "Selection from the reference file wasn't the same size as\n"
             << "the selection from the main system.  You must be able to use\n"
             << "the same selection string on both systems."
             << endl;
        exit(1);
        }
    sel.copyCoordinatesFrom(ref_sel);
    }
else if ( !(sel[0]->checkProperty(Atom::coordsbit)) )
    {
    traj->readFrame(0);
    traj->updateGroupCoords(system);
    }

bool use_periodicity_for_reference = topts->use_periodicity;
if (topts->use_periodicity && !system.isPeriodic())
    {
    use_periodicity_for_reference = false;
    cerr << "Warning: you requested periodicty, but the reference structure is not periodic" << endl;
    cerr << "Periodicity will _not_ be used when computing the reference contacts, " << endl;
    cerr << "but _will_ be used for the trajectory frames." << endl;
    }

// Compute the centers of mass of the selections
uint num_residues = residues.size();
vector<GCoord> centers_of_mass(num_residues);
for (uint i=0; i<num_residues; i++)
    {
    centers_of_mass[i] = residues[i].centerOfMass();
    }

vector<vector<uint> > contacts;
vector<uint>total_contacts_per_residue(num_residues);
vector<uint>contacts_per_residue(num_residues);

GCoord box = system.periodicBox();

int step = 1;
if (topts->exclude_consecutive)
    {
    step = 2;
    }
// Find contacts within the threshold distance
for (uint i=0; i<num_residues-step; i++)
    {
    for (uint j=i+step; j< num_residues; j++)
        {
        GCoord diff = centers_of_mass[j] - centers_of_mass[i];
        if (use_periodicity_for_reference)
            {
            diff.reimage(box);
            }
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

            // Store the total number of contacts for each residue
            if (topts->do_per_residue)
                {
                total_contacts_per_residue[i] += 1;
                total_contacts_per_residue[j] += 1;
                }
            }
        }
    }

// Number of native contacts, as a float because we'll need
// to do floating point arithmatic with it later anyway
float num_native_contacts = (float) contacts.size();
cout << "# Total native contacts: " << num_native_contacts << endl;

bool is_periodic;
if (topts->use_periodicity && traj->hasPeriodicBox())
    {
    is_periodic = true;
    }
else if (topts->use_periodicity && !(traj->hasPeriodicBox()))
    {
    cerr << "Warn: you requested periodicity, but your trajectory isn't periodic." << endl;
    cerr << "The calculation will proceed _ignoring_ periodicity." << endl;
    is_periodic = false;
    }

// Loop over structures in the trajectory
vector<vector<uint> >::iterator p;
int frame = 0;
while (traj->readFrame())
    {
    traj->updateGroupCoords(system);
    box = system.periodicBox();

    // Loop over contacts from the native structure
    int num_contacts = 0;
    for (p=contacts.begin(); p!= contacts.end(); p++)
        {
        uint r1 = p->at(0);
        uint r2 = p->at(1);
        GCoord c1 = residues[r1].centerOfMass();
        GCoord c2 = residues[r2].centerOfMass();
        GCoord diff = c2 - c1;
        if (is_periodic)
            {
            diff.reimage(box);
            }
        if (diff.length2() <= cut2)
            {
            num_contacts++;
            if (topts->do_output) output << "1\t";
            if (topts->do_per_residue)
                {
                contacts_per_residue[r1]++;
                contacts_per_residue[r2]++;
                }
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

// Output total contacts per residue
if (topts->do_per_residue)
    {
    ofstream per_residue_stream;
    per_residue_stream.open(topts->per_residue_filename.c_str());
    if (!per_residue_stream.is_open())
        {
        throw(runtime_error("couldn't open per-residue output file"));
        }
    per_residue_stream << "# Residue\tAveContacts\tTotalContacts" << endl;
    for (uint i=0; i < num_residues; ++i)
        {
        int real_residue = (residues[i][0])->resid();
        per_residue_stream << real_residue << "\t";
        if (total_contacts_per_residue[i])
            {
            double ave = static_cast<double>(contacts_per_residue[i])/
                (total_contacts_per_residue[i] * frame);
            per_residue_stream << ave;
            }
        else
            {
            per_residue_stream << "-1";
            }
        per_residue_stream << "\t" << total_contacts_per_residue[i] << endl;
        }
    per_residue_stream.close();
    }

}
