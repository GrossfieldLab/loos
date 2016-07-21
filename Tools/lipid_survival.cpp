/*
  Lipid lifetime
  
  (c) 2009 Joshua Horn, Grossfield Lab
      2016 Alan Grossdfield
  Department of Biochemistry and Biophysics
  University of Rochster School of Medicine and Dentistry

  Given a lipid in contact with a protein at time t, what is the prob
  that the lipid will be in contact at time t+dt?

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
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <numeric>

using namespace std;
using namespace loos;
namespace opts = loos::OptionsFramework;
namespace po = boost::program_options;

typedef vector<AtomicGroup> vGroup;
typedef vector<vector <AtomicGroup> > list_vGroup; 


string fullHelp(void) {
    string s =
"\n"
"SYNOPSIS\n"
"\n"
"    Compute the survival probability for a target molecule type around a probe\n"
"\n"
"DESCRIPTION\n"
"\n"
"This tool is used to calculate the survival probability for some kind of\n"
"probe molecule (e.g. a lipid) around a target molecule (e.g. a protein).\n"
"\n"
"The survival probability is the probability that, if the probe molecule \n"
"is \"bound\" at time t, it will also be bound at time t+delta t.  When\n"
"plotted as a function of delta t, this probability will decay from 1 to \n"
"0, and can generally be fit by a sum of exponentials.\n"
"\n"
"In general, one would more commonly use a correlation function here.\n"
"However, if the decay time is on the same timescale as your simulation,\n"
"the correlation function can go negative at long times (essentially saying\n"
"that lipids found at the protein surface early in the simulation are\n"
"unlikely to be present at the end, as opposed to being random).\n"
"Correlation functions with negative values are a pain to work with, so\n"
"we use survival probabilty as a convenient proxy.\n"
"\n"
"NOTE: The name \"survival probability\" could be slightly misleading; "
"      The quantity plotted is \n"
"      P_bound(t+dt|t)\n"
"      and does _not_ imply that the molecule was bound continuously \n"
"      during that interval.\n"
"\n"
"EXAMPLE\n"
"   lipid_lifetime --maxdt 2500 --probe 'segid == \"PROT\" && !hydrogen' --target 'resname == \"SDPE\" && name =~ \"C2\\d+\"' struct.pdb struct.dcd\n"
"\n"
"This will compute the correlation out to 2500 frames, looking for contacts\n"
"between heavy atoms in PROT and the saturated carbons in SDPE lipids.\n"
"Each lipid is considered separately, and the results are averaged over all\n"
"selected lipids.\n";
    return(s);
}

class ToolOptions : public opts::OptionsPackage {
public:

    ToolOptions() 
        {
        }
    void addGeneric(po::options_description& o) 
        {
        o.add_options()
      ("probe,p", po::value<string>(&protein_selection), "Main selection (e.g. protein)")
      ("target,t", po::value<string>(&lipid_selection), "Target selection (e.g. lipids)")
      ("cutoff,c", po::value<double>(&cutoff)->default_value(6.0), "Cutoff distance for contact")
      ("maxdt,m", po::value<uint>(&maxdt)->default_value(1000), "Maximum dt to compute")
      ("reimage,r", po::value<bool>(&reimage)->default_value(false), "Perform contact calculations considering periodicity")
      ;
        }

    string protein_selection;
    string lipid_selection;
    double cutoff;
    uint maxdt;
    bool reimage;
};

int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* basic = new opts::BasicOptions(fullHelp());
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory; 
  ToolOptions* topts = new ToolOptions;
  opts::AggregateOptions options;
  
  options.add(basic).add(tropts).add(topts);

  
  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;

  AtomicGroup protein = selectAtoms(model, topts->protein_selection);

  cout << "# " << hdr << endl;
  // selections for targets
  AtomicGroup lipid = selectAtoms(model, topts->lipid_selection);
  vGroup lipids = lipid.splitByMolecule();


  vector < vector <double> > contacts;
  // initialize vector of timeseries to store 
  for (unsigned int i = 0; i < lipids.size(); i++)
      {

        vector <double> tmp (traj->nframes(), 0.0);
        contacts.push_back(tmp);

      }

int frame_count = 0;
while (traj->readFrame()) 
    {
    traj->updateGroupCoords(model);
    GCoord box = model.periodicBox();
    
    for (uint j=0; j < contacts.size(); j++)
        {
        bool contact = false;
        if (topts->reimage) 
            {
            contact = lipids[j].contactWith(topts->cutoff, protein, box);
            }
        else
            {
            contact = lipids[j].contactWith(topts->cutoff, protein);
            }
    
        if (contact)
            {
            contacts[j][frame_count] = 1.0;
            }
        }
      frame_count++;
      }

/* Probability Calculations
 */

// loop over dt
//
cout << "0\t1.00" << endl;
for (unsigned int t = 1; t < topts->maxdt; t++)
    {

    double bound_tmp = 0.0;
    double total_tmp = 0.0;
    // loop over lipids
    for (unsigned int i = 0; i < contacts.size(); i++)
        {

        // loop through values for that lipid's contacts
        for (unsigned int j = 0; j < (contacts[i].size() - t); j++)
            {
            
            // if current value is 1
            if (contacts[i][j] == 1)
                {
                    if ((contacts[i][j+t]) == 1)
                        {

                        bound_tmp++;

                        }

                        total_tmp++;
                }

            }

        }

    double prob_tmp = bound_tmp/total_tmp;

    cout << t << "\t" << prob_tmp << endl;
    }
} 

