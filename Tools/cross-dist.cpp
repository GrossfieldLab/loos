/*
  cross-dist.cpp

  (c) 2012 Alan Grossfield 
           Department of Biochemistry
           University of Rochster School of Medicine and Dentistry

    Tool to compute the distribution of crossing angles between chains

*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011 Tod D. Romo
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



double cutoff;
unsigned int num_bins;
bool use_cosine;


// The following conditional prevents this class from appearing in the
// Doxygen documentation for LOOS:
//
// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:

  // Change these options to reflect what your tool needs
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("num_bins", po::value<unsigned int>(&num_bins)->default_value(20), "Number of histogram bins")
      ("cutoff", po::value<double>(&cutoff)->default_value(10.0), "Distance cutoff for neighboring chains")
      ("use-cosine", po::value<bool>(&use_cosine)->default_value(false), "Histogram the cosine instead of the angle")
      ;
  }

  // The print() function returns a string that describes what all the
  // options are set to (for logging purposes)
  string print() const {
    ostringstream oss;
    oss << boost::format("num_bins=%d, cutoff=%f ") % num_bins % cutoff;
    return(oss.str());
  }

};
// @endcond
// ----------------------------------------------------------------

string fullHelpMessage(void)
    {
    string s = 
"\n"
"SYNOPSIS\n"
"\n"
"Compute the probability distribution crossing angles for a set of chains\n"
"\n"
"DESCRIPTION\n"
"\n"
"The purpose of this tool is to compute the distribution of crossing \n"
"angles and torsions for a set of chains.  Specifically, it takes a selection of \n"
"atoms, splits them into in individual chains based on connectivity,\n"
"and at each time point computes their centroids and principal axes.\n"
"If a pair of chains centroids are within a threshold distance, it \n"
"computes the angle between their first principle axes and histograms\n"
"it.  The absolute value of the angle is used, because the principal axis\n"
"calculation doesn't determine sign (meaning for a chain lying along the\n"
"x-axis you could get (1,0,0) or (-1,0,0). \n"
"\n"
"It also computes the torsion angle between the two chains, by generating\n"
"an extra point for each chain by stepping away from the centroid along\n"
"the principal axis.  In this case, the angle is mapped into the range\n"
"0-90 degrees, again because the principal axis calculation doesn't \n"
"determine sign.  As a result, the column with the torsion values will\n" 
"will always be zeroes above 90 degrees.\n"
"\n"
"The model file must contain connectivity information.\n"
"\n"
"Command-line options:\n"
"    --num_bins      number of bins in the histogram, which goes \n"
"                    0-180 deg, default = 20\n"
"    --cutoff        distance below which two chains are considered \n"
"                    neighbors, default = 10 ang\n"
"    --use-cosine    histogram the cosine of the angle instead of the angle,\n"
"                    which changes the histogram bounds to -1:1.\n"
"\n"
"EXAMPLE\n"
"\n"
"cross-dist --selection 'name =~ \"^C\\d+$\" && resname =~\"PALM|OLEO\"' namd.psf trj_1.dcd\n"
"\n"
"This example selects the PALM and OLEO chain carbons from a POPC bilayer, \n"
"and uses the default bin number and cutoff.\n"
"\n"
"The output would look like:\n"
"# cross-dist '--selection' 'name =~ \"^C\\d+$\" && resname =~\"PALM|OLEO\"' 'namd.psf' 'trj_1.dcd' - alan (Mon Apr  2 12:57:16 2012) {/home/alan/projects/LOOS/trunk/Packages/User} [2.0.0 120402]\n"
"# Number of chains: 360\n"
"# Total points = 332402  332402\n"
"# Ang   Ang     Tors\n"
"4.5     0.0602554       0.171936\n"
"13.5    0.128146        0.160808\n"
"22.5    0.156873        0.14607\n"
"(and more lines like this)\n"
"\n"
"The two numbers in the \"Total points\" line are the number of angles and\n"
"torsions used; if these aren't the same, something very strange has \n"
"happened.\n"
    ;
    return(s);
    }


int main(int argc, char *argv[]) {
  
  string header = invocationHeader(argc, argv);
  cout << "# " << header << endl;
  
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());

  // This tool can operate on a subset of atoms.  The BasicSelection
  // object provides the "--selection" option.
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::TrajectoryWithFrameIndices* tropts = 
                                new opts::TrajectoryWithFrameIndices;

  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts);

  if (!options.parse(argc, argv))
    exit(-1);

  // Pull the model from the options object (it will include coordinates)
  AtomicGroup model = tropts->model;
  
  // Pull out the trajectory and the list of frames
  pTraj traj = tropts->trajectory;
  vector<uint> frame_indices = tropts->frameList();


  // Select chains
  AtomicGroup all_chains = selectAtoms(model, sopts->selection);

  // Verify we have connectivity, then split into individual chains
  if (! all_chains.hasBonds())
      {
      cerr << "The selection doesn't appear to have any bonds, and " 
           << endl
           << "this program requires connectivity information." 
           << endl
           << "You need to use a model file that has bond information, "
           << endl
           << "a PSF or a PDB with CONECT records."
           << endl;
      exit(-1);
      }
           
  vector<AtomicGroup> chains = all_chains.splitByMolecule();

  cout << "# Number of chains: " << chains.size() << endl;

  double cutoff2 = cutoff * cutoff;

  // Allocate space to store the histogram of angles
  vector<uint> ang_hist(num_bins);
  vector<uint> tors_hist(num_bins);
  double hist_min, hist_max;
  if (use_cosine)
      {
      hist_min = 0.0;
      hist_max = 1.0;
      }
  else
      {
      hist_min = 0.0;
      hist_max = 180.0;
      }

  double bin_size = (hist_max - hist_min) / num_bins;
  uint ang_total = 0;
  uint tors_total = 0;

  // Allocate space to hold the chain centroids and principle axes
  vector<GCoord> centers(chains.size());
  vector<GCoord> paxes(chains.size());
  vector<GCoord> points(chains.size());

  // Loop over the trajectory
  for (uint i=0; i<frame_indices.size(); ++i)
    {
    // Read in the next frame
    traj->readFrame(frame_indices[i]);

    // Update the coordinates ONLY for the subset of atoms we're
    // interested in...
    traj->updateGroupCoords(all_chains);
    GCoord box = all_chains.periodicBox();

    // Compute the centroid and principle axes for each chain 
    for (uint i = 0; i < chains.size(); i++)
        {
        centers[i] = chains[i].centroid();
        // principalAxes returns a vector of GCoord, the first one
        // is the most significant vector
        paxes[i] = chains[i].principalAxes()[0];

        // JH 07-2013
        // Make sign of paxes vector significant by matching end-to-end vector
        GCoord end_to_end = (chains[i].getAtom(chains[i].size() - 1))->coords() - (chains[i].getAtom(0))->coords();
        if (end_to_end * paxes[i] < 0.0)
          paxes[i] *= -1.0;
        // end JH 07-2013

        points[i] = centers[i] + paxes[i];
        }

    double ang = 0.0;
    for (uint i = 0; i < chains.size()-1; i++)
        {
        for (uint j =i+1; j < chains.size(); j++)
            {
            double dist2 = centers[i].distance2(centers[j], box);
            if (dist2 < cutoff2)
                {
                // Compute the angle
                // principle axes are already unit length
                double c = paxes[i] * paxes[j]; 
                // make sure roundoff didn't put cos out of range
                c = min(max(-1.0, c),1.0); 
                if (use_cosine)
                    {
                    ang = c;
                    }
                else
                    {
                    ang = fabs(acos(c) * 180.0/M_PI);
                    }
                if (!( (ang <= hist_min) || (ang >= hist_max) ))
                    { 
                    uint bin = (uint)((ang-hist_min)/bin_size);
                    ang_hist[bin]++;
                    ang_total++;
                    }

                // Compute the torsion
                double t = Math::torsion(points[i], centers[i], 
                                   centers[j], points[j]);
                // reflect back into the first quadrant, which is
                // the only meaningful range
                // Yes, I know this is the slow way, but it's idiot 
                // proof...

                // JH 07-2013
                // now torsion can be between 0 and 180
                // ang = acos( fabs(cos(t*M_PI/180.0)) ) * 180.0/M_PI;
                ang = fabs(t);
                if (use_cosine)
                    {
                    ang = cos(ang*M_PI/180.0);
                    }

                if (!( (ang <= hist_min) || (ang >= hist_max) ))
                    {
                    uint bin = (uint)((ang-hist_min)/bin_size);
                    tors_hist[bin]++;
                    tors_total++;
                    }
                }
            }
        }
    
    }

// Output the angle histogram
cout << "# Total points = " << ang_total << "  " << tors_total << endl;
cout << "# Ang\tAng\tTors" << endl;
for (uint i=0; i < num_bins; i++)
    {
    double ang = bin_size * (i+0.5) + hist_min;
    double ang_prob = (double)ang_hist[i] / ang_total;
    double tors_prob = (double)tors_hist[i] / tors_total;
    cout  << ang << "\t" 
          << ang_prob << "\t" 
          << tors_prob 
          << endl;
    }

}
