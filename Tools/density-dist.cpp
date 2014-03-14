/*
  density-dist.cpp

  Compute the charge/mass/electron density along the z dimension of a system

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


#include <ctype.h>

#include <loos.hpp>


using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

enum CalculationType { ELECTRON, CHARGE, MASS };


// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() :
    symmetrize(false),
    window(0),
    calc_type_desc("electron")
  { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("zsymmetry", po::value<bool>(&symmetrize)->default_value(symmetrize), "Symmetric with respect to Z")
      ("type", po::value<string>(&calc_type_desc)->default_value(calc_type_desc), "Calculation type (mass, charge, electron)")
      ("window", po::value<uint>(&window)->default_value(window), "Window size (in frames) for time series (0 = disabled)")
      ;
  }

  void addHidden(po::options_description& o) {
    o.add_options()
      ("selections", po::value< vector<string> >(&selections), "selections");
  }

  void addPositional(po::positional_options_description& pos) {
    pos.add("selections", -1);
  }

  bool postConditions(po::variables_map& map) {
    selections.insert(selections.begin(), "all");

    if (toupper(calc_type_desc[0]) == 'C')
      calc_type = CHARGE;
    else if (toupper(calc_type_desc[0]) == 'E')
      calc_type = ELECTRON;
    else if (toupper(calc_type_desc[0]) == 'M')
      calc_type = MASS;
    else {
      cerr << "Error- unknown calculation type '" << calc_type_desc << "' (should be either charge, mass or electron)\n";
      return(false);
    }
    return(true);
  }

  string help() const {
    return(string(" [selection [selection ...]]"));
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("zsymmetry=%d, type='%s', window=%d, selections='") 
      % symmetrize
      % calc_type_desc
      % window;
    for (vector<string>::const_iterator i = selections.begin(); i != selections.end(); ++i)
      oss << *i << "'" << (i == selections.end()-1 ? "" : ",'");
    return(oss.str());
  }

  bool symmetrize;
  uint window;
  string calc_type_desc;
  vector<string> selections;
  CalculationType calc_type;
};

// @endcond

string fullHelpMessage(void)
    {
    string s =
       "\n" 
       " SYNOPSIS\n"
       "\n"
       " Compute the electron, mass, or charge density for the system and\n"
       " its components along the z-axis.\n"
       "\n"
       " DESCRIPTION\n"
       "\n"
       " The purpose of this tool is to computed the distribution of a system\n"
       " along the z-axis.  This is most useful for membrane systems, where \n"
       " the data provided is analogous to that from x-ray or neutron \n"
       " scattering.  By default the program computes the total distribution,\n"
       " but if 1 or more selections are given on the command line, the \n"
       " individual distributions for those selections are output as well.\n"
       " In addition, the program can measure the time dependence of the\n"
       " distribution, and can automatically symmetrize the distribution\n"
       " around z=0.\n"
       "\n"
       " If the box size fluctuates (e.g. this is a constant pressure or \n"
       " constant tension run), then the variation of the area in the x-y\n"
       " plane is taken into account.  The output units are FOO/Ang^3, where\n"
       " FOO is either mass in AMU or charge/electron density in electrons.\n"
       "\n"
       " Options\n"
       " --type      Type of distribution (mass, electron, or charge).  If the\n"
       "             system file provides this information, it is used.  If not,\n"
       "             there's a warning message and reasonable guesses are \n"
       "             provided.  Mass and electron densities are comparable to\n"
       "             the results of neutron and x-ray scattering experiments,\n"
       "             while charge densities can be used to compute the \n"
       "             electrostatic potential profile (see below).\n"
       " --zsymmetry symmetrize the distribution with respect to z=0.  This \n"
       "             assumes the trajectory has already been recentered such\n"
       "             that the membrane center is at z=0 (if not, you can do \n"
       "             this with recenter-trj or merge-traj).\n"
       " --skip      Number of frames to discard from the beginning of the\n"
       "             trajectory\n"
       "\n"
       " Options for time-dependent output\n"
       "\n"
       " If you wish to track the change in the distribution over time, you\n"
       " can specify the following options:\n"
       "\n"
       " --window    Integer specifying how often to output running averages, \n"
       "             in frames.\n"
       " --prefix    Name for the output files for windowed averages.  E.g. \n"
       "             --prefix foo would give output files foo_1.dat, foo_2.dat,\n"
       "             etc.  If prefix contains a directory name, the program\n"
       "             does not check to ensure that the directory exists, and\n"
       "             will fail with an error message if it doesn't.\n"
       "\n"
       " EXAMPLE\n"
       "\n"
       " density-dist --type=charge -- namd.psf merged_1ns.dcd -38 38 76 'resname ==\"PEGL\"' 'resname == \"PGGL\"' 'segid == \"BULK\"'\n"
       "\n"
       " This command line computes a charge density along the membrane normal,\n"
       " running from -38 to 38 angstroms, with 1 angstrom bins.  In addition \n"
       " to computing the full charge distribution, the charge distribution of\n"
       " 3 components is also computed, corresponding to 2 difference lipid\n"
       " headgroups and water. \n"
       "\n"
       " Note: the \"--\" after the --type argument is necessary to tell the\n"
       "       code to stop processing arguments as if they were command-line\n"
       "       flags.  If you don't include it, it will read the -38 as the\n"
       "       flag -3 with a value 8, and will choke.  \n"
       "\n"
       " The first few lines of output from this command will look like:\n"
       " # density-dist '--type=charge' '--' 'namd.psf' 'merged_1ns.dcd' '-38' '38' '76' 'resname ==\"PEGL\"' 'resname == \"PGGL\"' 'segid == \"BULK\"' - alan (Thu Mar  8 11:21:24 2012) {/home/alan/projects/analysis_tools/scripts} [1.7.5 120308]\n"
       " # Z	AllAtoms Set(1)  Set(2)  Set(3) \n"
       " -37.5	0.000144657	0	0	8.99952e-05	\n"
       " -36.5	-4.88093e-05	0	0	-0.000136131	\n"
       " -35.5	1.51166e-06	0	0	-7.14851e-05	\n"
       " -34.5	-4.04959e-05	0	0	-0.00014739	\n"
       " -33.5	-0.000119295	0	-1.66837e-08	-0.000223778	\n"
       " -32.5	0.000201665	0	-4.33823e-07	8.88924e-05\n"
       " (with more lines following)\n"
       "\n"
       " The first column is the center of the histogram in z, the second\n"
       " is the distribution using all of the atoms, and the final three \n"
       " columns correspond to the distribution of the three selections\n"
       " specified on the command line.\n"
       "\n"
       " If you wish to use the charge density to compute the elecstrostatic\n"
       " potential along the membrane normal, you can combine the output\n"
       " from the above command with the tool potential_profile.py.  See \n"
       " the fullhelp message for that tool for more details.\n"
        ;
    return(s);
    }



int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);

  // Options handling...
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::OutputPrefix* popts = new opts::OutputPrefix;
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;

  // The required options could have been folded into ToolOptions, but
  // using RequiredArguments saves us from having to handle the help()
  // and print() methods...
  opts::RequiredArguments* ropts = new opts::RequiredArguments;
  ropts->addArgument("minz", "min-z");
  ropts->addArgument("maxz", "max-z");
  ropts->addArgument("nbins", "number-of-bins");
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(popts).add(tropts).add(ropts).add(topts);
  if (!options.parse(argc, argv)) {
    cerr << endl;
    cerr << "**Important note**\nYou must place '--' on the command line AFTER\n";
    cerr << "the options if you are going to use a negative Z argument, i.e.\n";
    cerr << "density-dist --type charge -- foo.pdb foo.dcd -40 40 40\n";
    exit(-1);
  }

  double min_z = parseStringAs<double>(ropts->value("minz"));
  double max_z = parseStringAs<double>(ropts->value("maxz"));
  uint nbins = parseStringAs<uint>(ropts->value("nbins"));
  vector<string> selections = topts->selections;

  AtomicGroup system = tropts->model;
  pTraj traj = tropts->trajectory;
  // End of options  
  cout << "# " << hdr << endl;

  // density from each selection
  vector<AtomicGroup> subsets;
  for (vector<string>::iterator i = selections.begin(); i != selections.end(); ++i)
    subsets.push_back(selectAtoms(system, *i));

  // Verify properties...
  for (vector<AtomicGroup>::iterator i = subsets.begin(); i != subsets.end(); ++i) {
    bool ok;
    switch(topts->calc_type) {
    case CHARGE: ok = i->allHaveProperty(Atom::chargebit); break;
    case ELECTRON: ok = (i->allHaveProperty(Atom::anumbit) && i->allHaveProperty(Atom::chargebit)); break;
    case MASS: ok = i->allHaveProperty(Atom::massbit); break;
    default: cerr << "Internal Error- unknown calculation type\n"; exit(-1);
    }
    if (!ok) {
      cerr << "***WARNING***\n";
      cerr << "The system is missing properties required for calculation type.\n";
      cerr << "Default values will be used where possible.\n";
      cerr << "This may result in incorrect or absurd values.\n";
      break;
    }
  }
  
  double bin_width = (max_z - min_z) / nbins;

  // Create and zero out the total charge_distribution
  vector< vector<double> > dists(subsets.size(), vector<double>(nbins, 0.0));
  // If we do windowed time series too, we'll need to zero out as we go, which 
  // means we'll need separate storage for the sum
  vector< vector<double> > cum_dists(subsets.size(), vector<double>(nbins, 0.0));

  // Note: the equillibration frames are already skipped by opts::BasicTrajectory
  
  // loop over the remaining frames
  uint frame = 0;
  while (traj->readFrame()) {
    // update coordinates
    traj->updateGroupCoords(system);

    // Compute the bin volume for normalization purposes
    GCoord box = system.periodicBox();
    double bin_volume = bin_width * box.x() * box.y();

    // Loop over the subsets and compute charge distribution.
    // (first set is all atoms)
    for (uint i=0; i<subsets.size(); i++) {
      double weight;
      for (AtomicGroup::iterator atom = subsets[i].begin(); 
                                 atom != subsets[i].end(); 
                                 ++atom) {
        switch(topts->calc_type) {
        case CHARGE: weight = (*atom)->charge(); break;
        case MASS: weight = (*atom)->mass(); break;
        case ELECTRON: weight = (*atom)->atomic_number() - (*atom)->charge(); break;
        default:
          cerr << "ERROR- unknown calculation type\n";
          exit(-1);
        }
        double z = (*atom)->coords().z();
        if (topts->symmetrize)
          z = abs(z);
        
        if ( (z > min_z) && (z < max_z) ) {
          uint bin = (int)( (z - min_z) / bin_width );
          dists[i][bin] += weight/bin_volume;
        }
      }
    }
    ++frame;

    // If windowed time series were requested, output them here
    if (topts->window && (frame % topts->window == 0)) {
        uint output_val = frame / topts->window;

        // build the output file name
        char piece[128];
        snprintf(piece, 128, "_%d.dat", output_val);
        string file_name = popts->prefix + string(piece);
        ofstream outfile(file_name.c_str());
        if (!outfile.is_open() ) {
            throw(runtime_error("couldn't open output file"));
            }

        // write the header
        outfile << "# Z\t";
        outfile << "AllAtoms\t";
        for (uint i=1; i<subsets.size(); i++) {
          outfile << " Set(" << i <<") "; 
        }
        outfile << endl;

        // Normalize by the number of frames and
        // output the average density
        for (uint i=0; i<nbins; i++) {
            double z = (i+0.5)*bin_width + min_z;
            outfile << z << "\t";
            for (unsigned int j=0; j<subsets.size(); j++) {
                double c = dists[j][i] / topts->window;
                outfile << c << "\t";
            }
            outfile << endl;
        }

        // Add the windowed density to the cumulative density and
        // zero out the distributions
        
        for (uint i=0; i<dists.size(); i++) {
            for (uint j=0; j<nbins; j++) {
                cum_dists[i][j] += dists[i][j];
            }

            // zero out the windowed dist
            dists[i].assign(dists[i].size(), 0.0);
        }
        
    }


  }

  // Normalize by the number of frames and output the average charge density
  cout << "# Z\tAllAtoms";
  for (uint i=1; i<subsets.size(); i++) {
    cout << " Set(" << i <<") "; 
  }
  cout << endl;

  // If we didn't output windowed time series, we never used cum_dists, so we need 
  // to copy it over here
  // If we did, we may have some data in dists that hasn't been added to cum_dists yet
  if (!topts->window) {
      cum_dists = dists;
  }
  else if (frame % topts->window != 0){
     for (uint i=0; i<dists.size(); i++) {
         for (uint j=0; j<nbins; j++) {
             cum_dists[i][j] += dists[i][j];
         }
     }
  }

  for (uint i=0; i<nbins; i++) {
    double z = (i+0.5)*bin_width + min_z;
    cout << z << "\t"; 
    for (uint j=0; j<subsets.size(); j++) {
      double c = cum_dists[j][i] / frame;
      
      cout << c << "\t";
    }
    cout << endl;
  }


  
}
