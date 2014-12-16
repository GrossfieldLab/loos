/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2013, Tod D. Romo
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


// @cond TOOL_INTERNAL
string fullHelpMessage(void) {
    string msg =
        "\n"
        "SYNOPSIS\n"
        "Convert a trajectory with N copies of a molecule into a longer one with only one copy\n"
        "\n"
        "DESCRIPTION\n"
        "\tGiven a trajectory that is L frames long and with N copies of a molecule (selection),\n"
        "this tool will write out a trajectory that is N*L frames long but with only one copy\n"
        "of the selection.  Each L-sized chunk of the trajectory contains a different copy\n"
        "of the same molecule.  For example, if a system has 100 frames and 4 identical peptides,\n"
        "then the output trajectory will be 400 frames long.  The first 100 frames will be\n"
        "the first peptide, the next 100 frames will be the second peptide, etc.  By default,\n"
        "the different molecules will not be centered.  The --center option will recenter\n"
        "the entire trajectory.\n"
        "\n"
        "\tserialize-selection can also be used to make a library of conformations from a\n"
        "trajectory.  The --pdbout option sets the prefix (using a printf-style format), and\n"
        "will cause the tool to write out a library of different PDB files rather than a\n"
        "DCD.  For membrane systems, it is convenient to canonicalize the orientation of the\n"
        "resulting library.  The --canon option will flip the selection depending on which\n"
        "leaflet it came frame (assuming the membrane normal points along the Z-axis).\n"
        "\n"
        "EXAMPLES\n"
        "\n"
        "\tserialize-selection --selection 'segid =~ \"PEP.\"' --prefix serial model.psf sim.dcd\n"
        "Extracts peptides with segid PEP0, PEP1, etc, into serial.pdb and serial.dcd\n"
        "\n"
        "\tserialize-selection --selection 'resname == \"POPC\"' --pdbout 'popc-%03d.pdb' --canon 1 model.psf sim.dcd\n"
        "Extracts POPC residues, centering them (by default, using all atoms in the molecule)\n"
        "and flipping them if they come from the lower leaflet.  The output is a library of PDB\n"
        "files named popc-000.pdb, popc-001.pdb, popc-002.pdb, etc.\n"
        "\n"
        "NOTES\n"
        "\tThe selection will be split into molecules either using connectivity (if present)\n"
        "or the segid.\n";

    return(msg);
}



class ToolOptions : public opts::OptionsPackage {
public:

    void addGeneric(po::options_description& o) {
        o.add_options()
	    ("renum", po::value<bool>(&renum)->default_value(true), "Renumber atomids in output PDB")
            ("pdbout", po::value<bool>(&pdb_output)->default_value(false), "Output is a library of PDBs (prefix must be a printf-style pattern)")
            ("center", po::value<string>(&center_selection)->default_value(""), "Selection to use for centering (empty selection does no centering)")
            ("canon", po::value<bool>(&canonicalize)->default_value(false), "Canonicalize orientation (for membrane peptides, flip orientation about X-axis if in lower leaflet.  Implies centering)");
    }


    bool postConditions(po::variables_map& vm) 
        {
            if (canonicalize && center_selection.empty()) {
                cerr << "Warning- canonicalization is turned on, but no centering selection provided.\n";
                cerr << "         Centering entire molecule by default.\n";
                
                center_selection = string("all");
            }
            

            return(true);
        }
    
    
    string print() const {
        ostringstream oss;

        oss << boost::format("pdbout=%d,center='%s',canon=%d") % pdb_output % center_selection % canonicalize;
        return(oss.str());
    }

    
    bool pdb_output;
    string center_selection;
    bool canonicalize;
    bool renum;
};





// Base clase for determining how the output is handled...
class Outputter 
{
public:
    Outputter(const string& prefix, const string& hdr) 
        : _prefix(prefix), _hdr(hdr) 
        {}
    

    virtual void writeFrame(const AtomicGroup& structure) =0;

    virtual ~Outputter() 
        {}
    
protected:
    string _prefix, _hdr;
};



/*
 * Writes output as a set of PDBs (one per frame).  The name of each
 * PDB is determined by the output prefix.  This object tracks the
 * number of frames written and uses this to generate the output name.
 */

class PDBOutput : public Outputter 
{
public:
    PDBOutput(const string& prefix, const string& hdr) 
        : Outputter(prefix, hdr), _count(0) 
        {

            // Check to see if the format string really generates unique names...
            ostringstream ossa, ossb;
            ossa << boost::format(prefix) % 1;
            ossb << boost::format(prefix) % 2;


            if (ossa.str() == ossb.str()) {
                cerr << "Error- output prefix needs to be a printf-style format string\n"
                    "       when using pdb output mode, e.g. 'foo%05d.pdb'.\n";
                exit(-10);
            }
            
        }
    
    void writeFrame(const AtomicGroup& structure) 
        {
            PDB pdb = PDB::fromAtomicGroup(structure);
            pdb.remarks().add(_hdr);

            ostringstream oss;
            oss << boost::format(_prefix) % _count;
            ofstream ofs(oss.str().c_str());
            if (ofs.fail()) {
                cerr << "Error- failed to open '" << oss.str() << "' for output.\n";
                exit(-11);
            }
            
            ofs << pdb;

            ++_count;
        }

private:
    uint _count;
};



/*
 * Writes the output as a Trajectory.  The first frame that is passed is used
 * to generate a PDB model.  The prefix name is used to name the
 * output files.
 */

class TrajOutput : public Outputter 
{
public:
    
    TrajOutput(const string& prefix, const string& type, const bool append, const string& hdr, const bool renum) 
        : Outputter(prefix, hdr), _first_frame(true), _renum(renum)
        {
            _traj = createOutputTrajectory(prefix + "." + type, append);
            _traj->setComments(hdr);
        }

    void writeFrame(const AtomicGroup& structure) 
        {
            if (_first_frame) {
		AtomicGroup output_model = structure.copy();
		if (_renum)
		    output_model.renumber();
		
                PDB pdb = PDB::fromAtomicGroup(output_model);
                pdb.remarks().add(_hdr);
                string pdbname = _prefix + ".pdb";
                ofstream ofs(pdbname.c_str());
                if (ofs.fail()) {
                    cerr << "Error- failed to open '" << pdbname << "' for output.\n";
                    exit(-11);
                }
                
                ofs << pdb;

                _first_frame = false;
            }
            _traj->writeFrame(structure);
        }

private:
    bool _first_frame;
    bool _renum;
    pTrajectoryWriter _traj;
};



// @endcond



int main(int argc, char *argv[]) 
{

    string hdr = invocationHeader(argc, argv);

    opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
    opts::OutputPrefix* popts = new opts::OutputPrefix;
    opts::BasicSelection* sopts = new opts::BasicSelection("all");
    opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
    opts::OutputTrajectoryTypeOptions* otopts = new opts::OutputTrajectoryTypeOptions;
    ToolOptions* topts = new ToolOptions;

    opts::AggregateOptions options;
    options.add(bopts).add(popts).add(sopts).add(tropts).add(otopts).add(topts);
    if (!options.parse(argc, argv))
        exit(-1);
    
    AtomicGroup subset = selectAtoms(tropts->model, sopts->selection);
            
    vector<AtomicGroup> molecules;
    if (tropts->model.hasBonds())
        molecules = subset.splitByMolecule();
    else
        molecules = subset.splitByUniqueSegid();

    // Simple safety check for bad selections/connectivity
    for (uint i=1; i<molecules.size(); ++i)
        if (molecules[i].size() != molecules[0].size()) {
            cerr << "Error- molecule #" << i << " has a different size than the first one.\n";
            cerr << "       Check your selection and try again.\n";
            exit(-10);
        }
    
    AtomicGroup outgroup = molecules[0].copy();

    // Figure out how to center
    bool do_center = !topts->center_selection.empty();
    AtomicGroup centering_subset = outgroup;
    if (do_center)
        centering_subset = selectAtoms(outgroup, topts->center_selection);
    
    // Set output type
    Outputter* output = 0;
    if (topts->pdb_output)
        output = new PDBOutput(popts->prefix, hdr);
    else
        output = new TrajOutput(popts->prefix, otopts->type, otopts->append, hdr, topts->renum);
    
    vector<uint> frames = tropts->frameList();


    for (uint k=0; k<molecules.size(); ++k) {
	for (uint j=0; j<frames.size(); ++j) {
	    tropts->trajectory->readFrame(frames[j]);
	    tropts->trajectory->updateGroupCoords(tropts->model);
	    
	    for (uint i=0; i<outgroup.size(); ++i)
		outgroup[i]->coords(molecules[k][i]->coords());
	    
            if (topts->canonicalize) {
                GCoord c = outgroup.centroid();
                if (c.z() < 0.0) {
                    outgroup.translate(-c);
                    outgroup.rotate(GCoord(1,0,0), 180.0);
                }
            }
            
            if (do_center) {
                GCoord c = centering_subset.centroid();
                outgroup.translate(-c);
            }
            
	    output->writeFrame(outgroup);
	}
    }
    
}
