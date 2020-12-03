/*
 *  Compute membrane property distribution about a protein
 *  Can compute chain molecular order parameter, tilt vector, or density
 *
 *  Alan Grossfield
 *  (adapted from some code by Josh Horn, adapted from some stuff I wrote)
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
#include <string>
#include "membrane_map.hpp"

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

enum CalcType { DENSITY, ORDER, HEIGHT, VECTOR };

class ToolOptions: public opts::OptionsPackage
{
public:
    void addGeneric(po::options_description& o)
        {
        o.add_options()
            ("xmin", po::value<double>(&xmin)->default_value(-50), "x histogram range")
            ("xmax",  po::value<double>(&xmax)->default_value(50), "x histogram range")
            ("xbins",  po::value<uint>(&xbins)->default_value(50), "x histogram bins")
            ("ymin", po::value<double>(&ymin)->default_value(-50), "y histogram range")
            ("ymax",  po::value<double>(&ymax)->default_value(50), "y histogram range")
            ("ybins",  po::value<uint>(&ybins)->default_value(50), "y histogram bins")
            ("calc", po::value<string>(&calc_type)->default_value(string("density")), "property to calculate (density, height, order, vector)")
            ("upper-only", "Map only the upper leaflet")
            ("lower-only", "Map only the lower leaflet")
            ("ref-structure", po::value<string>(&reference_filename), "Align to an external structure instead of the first frame")
            ("target-selection", po::value<string>(&target_selection), "Selection to use to calculate property")
            ("align-selection", po::value<string>(&align_selection), "Selection used to align the system")
            ;
        }

    bool postConditions(po::variables_map& vm)
        {
        if (calc_type.compare(string("density"))==0)
            {
            type = DENSITY;
            }
        else if (calc_type.compare(string("height"))==0)
            {
            type = HEIGHT;
            }
        else if (calc_type.compare(string("order"))==0)
            {
            type = ORDER;
            }
        else if (calc_type.compare(string("vector"))==0)
            {
            type = VECTOR;
            }
        else
            {
            cerr << "Error: unknown calculation type '" << calc_type
                 << "' (must be density, height, order, or vector)"
                 << endl;
            return(false);
            }

        upper_only = false;
        lower_only = false;
        if (vm.count("upper-only"))
            {
            upper_only = true;
            }
        if (vm.count("lower-only"))
            {
            lower_only = true;
            }
        if (upper_only && lower_only)
            {
            cerr << "Can't specify --upper-only and --lower-only at the same time"
                 << endl;
            exit(-1);
            }

        if (vm.count("align_selection"))
            {
            has_align = true;
            }
        else
            {
            has_align = false;
            }
        return(true);
        }

    double xmin, xmax, ymin, ymax;
    uint xbins, ybins;
    string calc_type;
    string reference_filename;
    string align_selection;
    string target_selection;
    CalcType type;
    bool upper_only;
    bool lower_only;
    bool has_align;
};


string fullHelpMessage(void)
    {
    string msg =
"\n"
"SYNOPSIS\n"
"\n"
"Compute the distribution of a membrane physical property about a membrane\n"
"protein.\n"
"\n"
"DESCRIPTION\n"
"The purpose of this tool is to compute one of a number of physical \n"
"properties on a 2D grid surrounding a membrane protein.  The properties \n"
"number density, height, molecular order parameter, and orientation vector,\n"
"although the code is written to make it easy to add other quantities (see \n"
"below).  The system is aligned against the coordinates of the selection \n"
"specified with --align-selection using the first unskipped frame as \n"
"reference, unless the --ref-structure option is given, in which case \n"
"the file specified there is used (--align-selection is still applied)\n"
"The alignment is performed in two dimensions, so that the \n"
"lipid bilayer is not tilted or shifted; it is assumed that the bilayer\n"
"normal is the z-axis, and that the bilayer center is at z=0\n"
"\n"
"By default, the target is treated as a single large entity.  The --splitby\n"
"flag will let you break it up into individual molecules, segments, or \n"
"residues.\n"
"\n"
"Options\n"
"--calc       The type of calculation to be performed.\n"
"             density: number density of the selection\n"
"             height: average z-position of the centroid of the selection\n"
"             order: molecular order parameter (see below)\n"
"             vector: orientation vector\n"
"\n"
"             The molecular order parameter is calculated using the \n"
"             principal axes of the selection; the 2nd and 3rd axes are\n"
"             averaged, and plugged into the standard 0.5 (3 cos^2 - 1) \n"
"             formula, relative the z-axis.  If the selection is a lipid\n"
"             chain, then the values are comparable to the ones seen for\n"
"             lipid order parameters.\n"
"\n"
"             The orientation vector is average of the xy components of the\n"
"             first principal axis of the selection.  Unlike the other \n"
"             options, which return scalars, this returns a 2D vector, \n"
"             which can be plotted in gnuplot using the \"with vector\" \n"
"             option.\n"
"\n"
"\n"
"\n"
"EXAMPLE\n"
"\n"
"membrane_map --xmin -30 --xmax 30 --ymin -30 --ymax 30 --xbins 30 --ybins 30 --splitby mol example.psf dark_ensemble_20.dcd --align-selection 'segid == \"RHOD\"' --target-selection 'resname == \"DHA\"'\n"
"\n"
"          This sets the histograms to run from -30:30 in x and y, with \n"
"          2 ang x 2 ang bins.  It uses the segment name RHOD to align the\n"
"          snapshots, and uses DHA chains as the targets.  Since no \n"
"          calculation type is specified, a number density is calculated. \n"
"          The DHA chains are split up on the basis of connectivity.\n"
"\n"
"If you wish to examine membrane properties in general (e.g. for a phase-\n"
"separated membrane with no protein) you can choose to not use an alignment\n"
"selection.  However, since domains may drift around during the simulation, \n"
"you may want to run the code on discrete ranges of frames rather than just \n"
"averaging over the whole trajectory. For example, you could modify the \n"
"previous example to be:\n"
"\n"
"membrane_map --range 200:299 --xmin -30 --xmax 30 --ymin -30 --ymax 30 --xbins 30 --ybins 30 --splitby mol example.psf dark_ensemble_20.dcd --target-selection 'resname == \"DHA\"'\n"
"\n"
"This calculation would not perform any alignment, and would skip the first\n"
"200 frames, use the next 100 frames, then skip the rest of the trajectory.\n"
"\n"
"POTENTIAL COMPLICATIONS\n"
"\n"
"The code will break if the alignment and target selections overlap, \n"
"because the 2D alignment works by setting the z-coordinates of the \n"
"alignment selection to 0. \n"
"\n"
"In regions where there's no data (e.g. inside the region occluded by \n"
"whatever you're aligning to), the code outputs a value of -999999999.\n"
"chosen as the ANSI standard insane value.  The exception is if you're\n"
"doing a density calculation, in which case there's no divide by zero\n"
"and sanity reigns everywhere.\n"
"\n"
"The options --upper-only and --lower-only let you calculate properties\n"
"using only the upper and lower leaflets respectively.  The check is done\n"
"for each frame, so these options will handle the case where a component\n"
"is capable of flipping between leaflets on the MD timescale.  However,\n"
"the implementation assumes that the membrane has been previously centered\n"
"at z=0.  For obvious reasons, you can't specify both --upper-only and \n"
"--lower-only at the same time.  Note: the stated number of matching target\n"
"molecules output at the beginning of the run does not take this restriction\n"
"into account.\n"
"\n"
"IMPLEMENTING NEW QUANTITIES\n"
"\n"
"Implementing new quantities is quite easy as long as they return either\n"
"a scalar (double) or a vector (GCoord).  You begin by editing \n"
"membrane_map.hpp -- all of the classes that do the work are subclasses of \n"
"CalcProperty, so you'll need to create a new class analogous to the ones\n"
"already there, e.g. CalcDensity or CalcMolOrder.  You'll need to supply\n"
"a constructor that calls the CalcProperty constructor and a method called\n"
"calc that will do the calculation.  The last step of the calc method must\n"
"be calling the incr function, which will add the value into the histogram.\n"
"If the quantity being calculated is a density, you may need to supply a \n"
"normalize method as well (see CalcDensity for an example).\n"
"\n"
"Once the class is written, you just have to hook it in so that the binary\n"
"knows about it.  To do this, edit membrane_map.cpp.  First, you'll need to\n"
"add the new type to the enum CalcType (near the top of membrane_map.cpp).\n"
"Second, edit the postConditions method in the ToolOptions class to add \n"
"your calculation type (it's just a string of if-else-ifs).  Third, \n"
"edit the switch statement in the main body of the code to add your \n"
"new method for calculating (there's only 1 switch).  Finally, \n"
"if you intend anyone else to use your method, edit the documentation \n"
"string in addGeneric (also in ToolOptions) and the fullHelpMessage \n"
"function.\n"
"\n"
        "\n";
    return(msg);
    }




int main(int argc, char *argv[])
    {

    string hdr = invocationHeader(argc, argv);
    cout << "# " << hdr << endl;

    opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
    opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
    opts::BasicSplitBy *sopts = new opts::BasicSplitBy;
    opts::RequiredArguments* ropts = new opts::RequiredArguments;
    ToolOptions* topts = new ToolOptions;

    opts::AggregateOptions options;
    options.add(bopts).add(tropts).add(ropts).add(topts).add(sopts);
    if (!options.parse(argc, argv))
        {
        exit(-1);
        }

    AtomicGroup system = tropts->model;
    pTraj traj = tropts->trajectory;
    vector<uint> frames = tropts->frameList();
    traj->readFrame(frames[0]);
    traj->updateGroupCoords(system);



    AtomicGroup align_to;
    AtomicGroup reference;
    if (topts->has_align)
        {
        align_to = selectAtoms(system, topts->align_selection);
        if ((topts->reference_filename).length() > 0)
            {
            AtomicGroup reference_system = createSystem(topts->reference_filename);
            reference = selectAtoms(reference_system, topts->align_selection);
            }
        else
            {
            reference = align_to.copy();
            }
        }

    AtomicGroup apply_to = selectAtoms(system, topts->target_selection);

    vector<AtomicGroup> targets = sopts->split(apply_to);
    cout << "# Found " << targets.size() << " matching molecules" << endl;

    // Set up storage for our property.
    double xmin = topts->xmin;
    double xmax = topts->xmax;
    double ymin = topts->ymin;
    double ymax = topts->ymax;
    uint xbins = topts->xbins;
    uint ybins = topts->ybins;
    double xwidth = (xmax - xmin)/xbins;
    double ywidth = (ymax - ymin)/ybins;

    CalcPropertyBase *calculator;
    switch (topts->type)
        {
        case DENSITY:
            calculator = new CalcDensity(xbins, ybins, xwidth, ywidth);
            break;
        case ORDER:
            calculator = new CalcMolOrder(xbins, ybins);
            break;
        case HEIGHT:
            calculator = new CalcHeight(xbins, ybins);
            break;
        case VECTOR:
            calculator = new CalcOrientVector(xbins, ybins);
            break;
        default: // this can't happen, set in option handling
            cerr << "ERROR: unknown calculation type" << endl;
            exit(-1);
        }

    // We don't want the transformation to tilt the membrane, so we'll
    // zero out the z coordinates before using the alignment.
    // We'll do the same
    for (AtomicGroup::iterator i = reference.begin();
                               i!= reference.end();
                               ++i)
        {
        (*i)->coords().z() = 0.0;
        }

    // loop over frames in the trajectory
    for (uint i=0; i<frames.size(); ++i)
        {
        traj->readFrame(frames[i]);
        traj->updateGroupCoords(system);

        if (topts->has_align)
            {
            // zero out the alignment selections z-coordinate
            AtomicGroup align_to_flattened = align_to.copy();
            for (AtomicGroup::iterator j = align_to_flattened.begin();
                                       j!= align_to_flattened.end();
                                       ++j)
                {
                (*j)->coords().z() = 0.0;
                }


            // get the alignment matrix
            GMatrix M = align_to_flattened.superposition(reference);
            M(2,2) = 1.0;    // Fix a problem caused by zapping the z-coords...
            XForm W(M);

            // align the stuff we're goign to do the calculation on
            apply_to.applyTransform(W);
            }

        // Calculate something
        uint xbin, ybin;
        for (vector<AtomicGroup>::iterator j  = targets.begin();
                                           j != targets.end();
                                           ++j)
            {
            GCoord centroid = j->centroid();
            // Skip molecules outside the xy range of interest
            if ( (centroid.x() < xmin) || (centroid.x() > xmax) ||
                 (centroid.y() < ymin) || (centroid.y() > ymax)
               )
                {
                continue;
                }
            // If the user chose to look at only one leaflet,
            // skip molecules in the opposite leaflet.
            // Note: this assumes that the membrane is centered at z=0
            else if ((centroid.z() > 0) && topts->lower_only)
                {
                continue;
                }
            else if ((centroid.z() < 0) && topts->upper_only)
                {
                continue;
                }
            else
                {
                xbin = (uint)(centroid.x() - xmin) / xwidth;
                ybin = (uint)(centroid.y() - ymin) / ywidth;
                }

            // do the work
            calculator->calc(*j, xbin, ybin);
            }

        }

    // loop over bins and dump out the values
    calculator->normalize(frames.size());

    cout << "# X\tY\tValue(s)" << endl;
    for (uint i = 0; i < xbins; ++i)
        {
        double xval = xmin + xwidth*i;
        for (uint j = 0; j < ybins; ++j)
            {
            double yval = ymin + ywidth*j;

            cout << xval << "\t"
                 << yval << "\t"
                 << calculator->print(i,j)
                 << endl;
            }
        cout << endl;
        }
    }
