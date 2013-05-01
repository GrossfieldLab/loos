// (c) 2013 Tod D. Romo, Grossfield Lab, URMC


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
        "\n"
        "\n"
        "DESCRIPTION\n"
        "\n"
        "EXAMPLES\n"
        "\n"
        "\n"
        "\n"
        "\n"
        "SEE ALSO\n"
        "\t\n\n";

    return(msg);
}



class ToolOptions : public opts::OptionsPackage {
public:

    void addGeneric(po::options_description& o) {
        o.add_options()
            ("pdbout", po::value<bool>(&pdb_output)->default_value(false), "Output is a library of PDBs (prefix must be a printf-style pattern)")
            ("center", po::value<bool>(&center_selection)->default_value(false), "Center the selection");
        
    }

    string print() const {
        ostringstream oss;

        oss << boost::format("pdbout=%d,center=%d") % pdb_output % center_selection;
        return(oss.str());
    }

    
    bool pdb_output, center_selection;
    
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
 * Writes the output as a DCD.  The first frame that is passed is used
 * to generate a PDB model.  The prefix name is used to name the
 * output files, i.e. prefix.pdb and prefix.dcd
 */

class DCDOutput : public Outputter 
{
public:
    
    DCDOutput(const string& prefix, const string& hdr) 
        : Outputter(prefix, hdr), _first_frame(true), _dcd(prefix + ".dcd") 
        {}

    void writeFrame(const AtomicGroup& structure) 
        {
            if (_first_frame) {
                PDB pdb = PDB::fromAtomicGroup(structure);
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
            _dcd.writeFrame(structure);
        }
private:
    bool _first_frame;
    DCDWriter _dcd;
};



// @endcond



int main(int argc, char *argv[]) 
{

    string hdr = invocationHeader(argc, argv);

    opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
    opts::OutputPrefix* popts = new opts::OutputPrefix;
        opts::BasicSelection* sopts = new opts::BasicSelection("all");
    opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
    ToolOptions* topts = new ToolOptions;

    opts::AggregateOptions options;
    options.add(bopts).add(popts).add(sopts).add(tropts).add(topts);
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

    // Set output type
    Outputter* output = 0;
    if (topts->pdb_output)
        output = new PDBOutput(popts->prefix, hdr);
    else
        output = new DCDOutput(popts->prefix, hdr);
    
    vector<uint> frames = tropts->frameList();
        
    for (uint k=0; k<frames.size(); ++k) {
        tropts->trajectory->readFrame(frames[k]);
        tropts->trajectory->updateGroupCoords(tropts->model);

        for (uint j=0; j<molecules.size(); ++j) {
            for (uint i=0; i<outgroup.size(); ++i)
                outgroup[i]->coords(molecules[j][i]->coords());

            if (topts->center_selection)
                outgroup.centerAtOrigin();
            
            output->writeFrame(outgroup);
        }
    }
}
