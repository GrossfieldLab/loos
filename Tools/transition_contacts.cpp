/*
  transition_contacts: compute the fraction of native contacts in a trajectory,
  based on an initial structure.
 
  This file is part of LOOS.
  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2013, Nick Leioats, Alan Grossfield
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

class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() :
    selection("!(segid == 'BULK' || segid == 'SOLV' || hydrogen)")
  { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("selection", po::value<string>(&selection)->default_value(selection), "Selection for calculation")
      ("cutoff", po::value<double>(&cutoff)->default_value(8.0), "Cutoff to use for defining contacts")
      ("source-selection", po::value<string>(&source_sel)->default_value(""), "Selection specific to source model")
      ("sink-selection", po::value<string>(&sink_sel)->default_value(""), "Selection specific to sink model")
      ("timeseries", po::value<string>(&timeseries)->default_value(""), "Report contacts as a timeseries")
      ("include-heavy", po::value<bool>(&leave_heavy)->default_value(false), "Include backbone and hydrogen atoms");
  }

  void addHidden(po::options_description& o) {
    o.add_options()
      ("source-model", po::value<string>(&source_model), "Source model")
      ("sink-model", po::value<string>(&sink_model), "Sink model");
  }

  void addPositional(po::positional_options_description& p) {
    p.add("source-model", 1);
    p.add("sink-model", 1);
  }

  string help() const { return("source-model sink-model"); }
  string print() const {
    ostringstream oss;
    oss << boost::format("sink-sel='%s', source-sel='%s', timeseries=%s")
      % cutoff
      % sink_sel
      % source_sel
      % timeseries
      % leave_heavy;
    return(oss.str());
  }

  bool postConditions(po::variables_map& map) {
    if (sink_sel.empty()) {
      sink_sel = selection; 
      cerr << "Warning: Using --selection for sink\n";
    }
    if (source_sel.empty()) {
      source_sel = selection; 
      cerr << "Warning: Using --selection for source\n";
    }
    
    return(true);
  }
  
  double cutoff;
  string sink_model, source_model;
  string selection;
  string sink_sel, source_sel;
  string timeseries;
  bool leave_heavy;
};

// @endcond


string fullHelpMessage(void) {
  string msg = 
    "\n"
    "SYNOPSIS\n"
    "\n"
    "Calculate the normalized transition between two structures.\n"
    "\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "This tool will calculate the number of unique contacts in\n"
    "a trajectory in order to asses its transition from a source\n"
    "structure to a sink structure.  Contacts are defined using\n"
    "two input structures.  The set of contacts which differ between\n"
    "the two input structures are then used in the calculation.\n"
    "This tool explicitly uses the side-chain centriods for the\n"
    "calculation. (Note: In this context the word \'unique\' is used\n"
    "to describe those contacts present in only the source OR the\n"
    "sink structure.)\n"
    "\n"
    "The result is a normalized count of the number of contacts\n"
    "broken and formed for each frame of the trajectory.\n"
    "The output is 4 columns:\n"
    "\tFrame           - trajectory frame number\n"
    "\tContacts broken - Number of unique source contacts broken\n"
    "\tContacts formed - Number of unique sink contacts formed\n"
    "\tTransition      - Normalized sum of broken and formed\n"
    "\n"
    "Optionally, a timeseries of each contact's state may be\n"
    "written out.  In this case 1's are written for formed\n"
    "contacts and 0's for broken contacts.  Additionally, a\n"
    "list of the contacts is output for reference.\n"
    "\n"
    "\n"
    "\n"
    //
    "EXAMPLE\n"
    "transition_contacts --cutoff 8 --selection \'segid==\"PROT\"\' model.pdb traj.dcd source.pdb sink.pdb\n"
    "\tHere we are calculating how far our simulation has\n"
    "\ttransitioned away from the structure in source.pdb\n"
    "\ttowards the structure in sink.pdb.   The selection\n"
    "\tspecifies the entire segid \"PROT\" is used in the\n"
    "\tcalculation.  In this example the same selection is\n"
    "\tapplied to both the source and sink models as well.\n"
    "\tAn 8 angstrom cutoff is used to define connectivity.\n"
    "\t\n"
    "\t\n"
    "transition_contacts --sink_sel 'resid==\"PROT\"' --source_sel 'resid==\"PROT\"' --cutoff 8 --selection \'segid==\"PROT\"\' model.pdb traj.dcd source.pdb sink.pdb\n"
    "\tSame as the above command, but now the selections\n"
    "\tfor the source and sink models are separately specified.\n"
    "\tIMPORTANT: It is your responsibility to ensure that\n"
    "\tthe atoms selected in all three models match.  This\n"
    "\ttool will run regardless, pushing residues onto a \n"
    "\tvector.  \n"
    "\t\n"
    "\t\n"
    "transition_contacts --timeseries output-time  --selection \'segid==\"PROT\"\' model.pdb traj.dcd source.pdb sink.pdb\n"
    "\tSame options as the first example, but now we\n"
    "\toutput the timeseries of contacts to the file\n"
    "\toutput-time in addition to the standard output.\n"
    "\t\n"
    "\t\n"
    "transition_contacts --include-heavy 1 --timeseries output-time  --selection \'segid==\"PROT\"\' model.pdb traj.dcd source.pdb sink.pdb\n"
    "\tSame as the example above, but now the calculation\n"
    "\tincludes backbone atoms and hydrogen atoms.  Where\n"
    "\tin previous examples these were excluded from the\n"
    "\tcalculation.\n"
    "\t\n"
    "\t\n"
    "SEE ALSO\n"
    "native_contacts -\n"
    "This tool calculates the changes in contacts when a 2nd structure\n"
    "is not available.\n"
    "\t\n"
    "\t\n"
    "\n";
  return(msg);
}


//### MAIN
int main (int argc, char *argv[]){
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);
  cout << "# " << hdr << endl;
  
  AtomicGroup model = tropts->model;
  AtomicGroup system = selectAtoms(model, topts->selection);
  AtomicGroup init_model = createSystem(topts->source_model);
  AtomicGroup init = selectAtoms(init_model, topts->source_sel);
  AtomicGroup final_model = createSystem(topts->sink_model);
  AtomicGroup final = selectAtoms(final_model, topts->sink_sel);

  double cutoff = topts->cutoff;
  pTraj traj = tropts->trajectory;

  // Check to see if user requested to leave these atoms...
  if (!topts->leave_heavy){
    system.remove(system.select(BackboneSelector()));
    system.remove(system.select(HydrogenSelector()));
    init.remove(init.select(BackboneSelector()));
    init.remove(init.select(HydrogenSelector()));
    final.remove(final.select(BackboneSelector()));
    final.remove(final.select(HydrogenSelector()));
  }else{
    cerr << "WARNING: Leaving backbone and hydrogen atoms!\n";
  }
  // Split into residues
  vector<AtomicGroup> start_residues = init.splitByResidue();
  vector<AtomicGroup> final_residues = final.splitByResidue();
  vector<AtomicGroup> residues = system.splitByResidue();

  if ((start_residues.size() != final_residues.size()) || (residues.size() != final_residues.size())) {
    cerr << boost::format("Error: The trajectory has %d residues, the source has %d residues, and sink has %d residues.\n") 
      % residues.size() 
      % start_residues.size() 
      % final_residues.size();
    cerr << "\tThe source and sink selections must have the same number of residues.\n";
    exit(-1);
  }
  
  double cut2 = cutoff*cutoff;
  string timeseries_outfile = topts->timeseries;
  ofstream ofs;
  if (!timeseries_outfile.empty()){
    ofs.open(timeseries_outfile.c_str());
    if (ofs.fail()) {
      cerr << "Error- failed to open '" << timeseries_outfile << "' for output.\n";
      exit(-11);
    }
    ofs << "# " << hdr << endl;
    ofs << "# Changed contacts list:\n";
    ofs << "# ------------------------------\n";
  }


  //######################################################################
  //### Build the master list of unique contacts
  //######################################################################
  vector<AtomicGroup> reslist;
  vector<pair<uint, uint> > formed_connection_list;
  vector<pair<uint, uint> > broken_connection_list;

  for (uint j = 0; j < residues.size()-1; ++j) {
    uint jk;
    bool any_connection_flag = false;
    GCoord cj = start_residues[j].centerOfMass();
    for (uint i = j+1; i < residues.size(); ++i) {
      GCoord ci = start_residues[i].centerOfMass();
      GCoord start_diff = cj - ci;

      // If formed in starting structure...
      if (start_diff.length2() <= cut2){
        GCoord fci = final_residues[i].centerOfMass();
        GCoord fcj = final_residues[j].centerOfMass();
        GCoord final_diff = fcj - fci;
        // ...and if broken in ending structure
        if (final_diff.length2() > cut2){
          // Check if residue[j]'s first connection
          // (we only want [j] listed once)
          if (!any_connection_flag){
            jk = reslist.size();
            reslist.push_back(residues[j]);
            any_connection_flag = true;
          }
          uint ik = reslist.size();
          reslist.push_back(residues[i]);

          pair<uint,uint> conn_ptr(jk, ik);
          // Add these to the broken pair list
          broken_connection_list.push_back(conn_ptr);
          // Document the residues in the broken list
          ofs << "# broken: " << final_residues[j].getAtom(0)->resid() << " "   
              << final_residues[i].getAtom(0)->resid() << endl;
        }
      }

      // If broken in starting structure...
      if (start_diff.length2() > cut2){
        GCoord fci = final_residues[i].centerOfMass();
        GCoord fcj = final_residues[j].centerOfMass();
        GCoord final_diff = fcj - fci;
        // ...and if formed in ending structure
        if (final_diff.length2() <= cut2){
          if (!any_connection_flag){
            jk = reslist.size();
            reslist.push_back(residues[j]);
            any_connection_flag = true;
          }
          uint ik = reslist.size();
          reslist.push_back(residues[i]);
          
          pair<uint,uint> conn_ptr(jk, ik);
          // Add to the formed pair master list
          formed_connection_list.push_back(conn_ptr);
          // Document the residues in the formed list
          ofs << "# formed: " << final_residues[j].getAtom(0)->resid() << " "   
              << final_residues[i].getAtom(0)->resid() << endl;
        }
      }
    }
  }


  //######################################################################
  //### Set up params for the calculation
  //######################################################################
  // Tally of the number of contacts for normalization  
  uint total_formed_contacts = formed_connection_list.size();
  uint total_broken_contacts = broken_connection_list.size();
  
  cout << "# Total differences: \n";
  cout << "# Contacts broken: " << total_broken_contacts << endl;
  cout << "# Contacts formed: " << total_formed_contacts << endl;
  cout << "# Frame \t broken \t formed \t total \n";
  cout << "#-------------------------------------------\n";


  //######################################################################
  //### Iterate over trajectory frames 
  //### Calc the number of unique contacts broken/formed
  //######################################################################
  int frame = tropts->skip;
  while (traj->readFrame()){
    traj->updateGroupCoords(system);
    uint number_broken = 0;
    uint number_formed = 0;
    
    // Build a GCoord vector of CoM's 
    // with indices corresponding to reslist
    vector<GCoord> centers;
    for (uint gi = 0; gi < reslist.size(); ++gi){
      GCoord temp_com = reslist[gi].centerOfMass();
      centers.push_back(temp_com);
    }

    // Do the actual calculation
    // For broken list
    for (uint i = 0; i < broken_connection_list.size(); ++i){

      GCoord first_current  = centers[broken_connection_list[i].first];
      GCoord second_current = centers[broken_connection_list[i].second]; 
      bool broken = first_current.distance2(second_current) > cut2;
      if (broken)
          ++number_broken;
      if (!timeseries_outfile.empty())
          ofs << (broken ? '0' : '1') << '\t';
    }

    // For formed list
    for (uint i = 0; i < formed_connection_list.size(); ++i){
      GCoord first_current  = centers[formed_connection_list[i].first];
      GCoord second_current = centers[formed_connection_list[i].second]; 
      bool formed = first_current.distance2(second_current) < cut2;

      if (formed)
          ++number_formed;
      if (!timeseries_outfile.empty())
          ofs << (formed ? '1' : '0') << '\t';
    }

    if (!timeseries_outfile.empty())
      ofs << endl;

  
    // Format and output the data
    float num_transition = (static_cast<float>(number_broken) + number_formed) 
        / (total_broken_contacts + total_formed_contacts);
    float percent_broken = static_cast<float>(number_broken) / total_broken_contacts;
    float percent_formed = static_cast<float>(number_formed) / total_formed_contacts;
    cout << frame << "\t" << percent_broken << "\t" << percent_formed << "\t" << num_transition << endl;
    frame++;
  }
}
