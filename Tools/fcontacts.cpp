/*
  fcontacts

  Computes the number of contacts between a probe group and a set of target groups...
*/


/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008-2013 Tod D. Romo
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
#include <limits>

using namespace std;
using namespace loos;
namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

typedef vector<AtomicGroup> vGroup;




// @cond TOOL_INTERNAL

string fullHelpMessage() {
    string s = 
        "\n"
        "SYNOPSIS\n"
        "Determine the number of contacts between a probe selection and multiple targets\n"
        "\n"
        "DESCRIPTION\n"
        "\tcontact-time can be used to find the number of putative contacts between\n"
        "a probe set of atoms and a number of different target sets of atoms.  contact-time\n"
        "counts the number of times a target atom is within a given shell about any\n"
        "probe atom.  A matrix is constructed where each target is a column and each row\n"
        "represents a time point in the trajectory.\n"
        "\n"
        "\tThe matrix  can be normalized in two ways: row or column.\n"
        "Row normalization gives the percentage contact between the probe\n"
        "and each target relative to all contacts.  Column normalization\n"
        "gives the percentage contact between the probe and each target\n"
        "relative to the maximum number of contacts against the respective\n"
        "target.\n"
        "\n"
        "\tThe autoself option splits the probe selection into a set of\n"
        "molecules based on segid.  It then computes the contacts between\n"
        "all of these molecules (excluding self-to-self) and includes this\n"
        "as an extra column in the output.  As an example, suppose\n"
        "you have a number of AMLPs in a membrane, each with a different\n"
        "segid (i.e. PE1, PE2, ...) and you want to find the percentage\n"
        "contacts between the AMLPs and PEGL, PGGL, and each other.  The\n"
        "command for this would be:\n"
        "\n"
        "contact-time --autoself=1 model.pdb traj.dcd  'segid =~ \"PE\\d+\"'\\\n"
        "      'resname == \"PEGL\"' and 'resname == \"PGGL\"'\n"
        "\n"
        "This will automatically generate a new set of targets based\n"
        "on the probe selection, splitting them into separate molecules\n"
        "based on their segid.  It then computes the unique pair-wise\n"
        "contacts between each AMLP.  The total number of self-contacts\n"
        "is then included as an extra column in the output.\n"
        "\n"
        "EXAMPLES\n"
        "\n"
        "\tcontact-time --inner 0 --outer 4.5 model.psf traj.dcd 'segid == \"PEPT\"'\\\n"
        "\t  'resname == \"PEGL\"' 'segid == \"BULK\"'\n"
        "This example counts the number of contacts within 4.5 angstroms of any\n"
        "PEGL atom with PEPT atoms, and any BULK atom with PEPT atoms.  Row\n"
        "normalization is used, so each row represents the percent contact of\n"
        "each target, e.g.. 20% PEGL and 50% BULK at time 10ns\n"
        "\n"
        "\tcontact-time --inner 0 --outer 4.5 --rownorm 0 --colnorm 1 \\\n"
        "\t  model.psf traj.dcd 'segid == \"PEPT\"' 'resname == \"PEGL\"' \\\n"
        "\t  'segid == \"BULK\"'\n"
        "This example is as above, but the matrix is normalized down a column.\n"
        "Here, the data would show that at time 10 ns, PEPT makes a 20% contact\n"
        "with PEGL (relative to the maximum contact with PEGL), and likewise\n"
        "for BULK\n"
        "\n"
        "NOTES\n"
        "\tBy default, contact-time uses a distance filter to eliminate\n"
        "target atoms that are too far to be considered when looking\n"
        "at each probe atom.  The padding for the radius used to\n"
        "exclude target atoms can be adjusted with the '--fastpad' option.\n"
        "In the unlikely event the filter causes problems, it can\n"
        "be disabled with '--fast=0'.\n";
  
    return(s);
}



class ToolOptions : public opts::OptionsPackage
{
public:
    ToolOptions() :
        inner_cutoff(0.0),
        outer_cutoff(4.0),
        pad(1.0),
        probe_selection(""),
        symmetry(true),
        auto_split(true)
        { }


    void addGeneric(po::options_description& o) {
        o.add_options()
            ("inner", po::value<double>(&inner_cutoff)->default_value(inner_cutoff), "Inner cutoff (ignore atoms closer than this)")
            ("outer", po::value<double>(&outer_cutoff)->default_value(outer_cutoff), "Outer cutoff (ignore atoms further away than this)")
            ("reimage", po::value<bool>(&symmetry)->default_value(symmetry), "Consider symmetry when computing distances")
            ("split", po::value<bool>(&auto_split)->default_value(auto_split), "Automatically split probe selection")
            ("pad", po::value<double>(&pad)->default_value(pad), "Padding for filtering nearby atoms");
    }

    void addHidden(po::options_description& o) {
        o.add_options()
            ("system", po::value<string>(&system_selection), "System selection")
            ("probe", po::value<string>(&probe_selection), "Probe selection")
            ("target", po::value< vector<string> >(&target_selections), "Target selections");
    }

    void addPositional(po::positional_options_description& p) {
        p.add("system", 1);
        p.add("probe", 1);
        p.add("target", -1);
    }

    bool check(po::variables_map& map) {
        if (target_selections.empty() || probe_selection.empty() || system_selection.empty())
            return(true);

        return(false);
    }


    string help() const {
        return("system probe target [target ...]");
    }

    string print() const {
        ostringstream oss;
        oss << boost::format("inner=%f,outer=%f,reimage=%d,autosplit=%d,pad=%f,system='%s',probe='%s',targets=")
            % inner_cutoff
            % outer_cutoff
            % symmetry
            % auto_split
            % pad
            % system_selection
            % probe_selection;

        for (uint i=0; i<target_selections.size(); ++i)
            oss << "'" << target_selections[i] << "'" << (i == target_selections.size() - 1 ? "" : ",");

        return(oss.str());
    }


    double inner_cutoff, outer_cutoff, pad;
    string probe_selection, system_selection;
    bool symmetry, auto_split;
    vector<string> target_selections;
};



// Check for atom equality only through atomid...

struct IdEquals : public binary_function<pAtom, pAtom, bool> 
{
    bool operator()(const pAtom& a, const pAtom& b) const 
        {
            return(a->id() == b->id());
        }
};


typedef vector< vector<double> >      FContactsList;


// @endcond



// Return a list of target atoms that are in contact with probe
AtomicGroup contacts(const AtomicGroup& probe, const AtomicGroup& target, const double inner_radius, const double outer_radius, const bool symmetry) {
  
    double or2 = outer_radius * outer_radius;
    double ir2 = inner_radius * inner_radius;

    GCoord box = target.periodicBox();
    AtomicGroup contacting_atoms;
    

    for (AtomicGroup::const_iterator j = probe.begin(); j != probe.end(); ++j) {
        GCoord v = (*j)->coords();
    
        for (AtomicGroup::const_iterator i = target.begin(); i != target.end(); ++i) {
            GCoord u = (*i)->coords();
            double d = symmetry ? v.distance2(u, box) : v.distance2(u);
            if (d >= ir2 && d <= or2)
                contacting_atoms.append(*i);
        }
    }

    return(contacting_atoms);
}


AtomicGroup pickNearbyAtoms(const AtomicGroup& probe, const AtomicGroup& target, const double radius, const bool symmetry) {

    GCoord c = probe.centroid();
    GCoord box = probe.periodicBox();
    double maxrad = probe.radius() + radius;
    maxrad *= maxrad;

  
    AtomicGroup nearby;
    nearby.periodicBox(probe.periodicBox());
    for (AtomicGroup::const_iterator i = target.begin(); i != target.end(); ++i) {
        double d = symmetry ? c.distance2((*i)->coords(), box) : c.distance2((*i)->coords());
        if (d <= maxrad)
            nearby.append(*i);
    }

    return(nearby);
}



vector<double> fractionContactsToProbe(const AtomicGroup& probe,
                                       const AtomicGroup& nearby,
                                       const vGroup& targets,
                                       const double inner_radius,
                                       const double outer_radius,
                                       const bool symmetry)
{

    // First, find which nearby atoms are actually in contact...
    AtomicGroup nearby_contacts = contacts(probe, nearby, inner_radius, outer_radius, symmetry);
    AtomicGroup myself = nearby_contacts.intersect(probe, IdEquals());
    nearby_contacts.remove(myself);
    
    vector<double> fracts(targets.size(), 0.0);
    for (uint i=0; i<targets.size(); ++i) {
        AtomicGroup target_nearby = nearby_contacts.intersect(targets[i], IdEquals());
        fracts[i] = static_cast<double>(target_nearby.size()) / nearby_contacts.size();
    }
    return(fracts);
}

        
FContactsList fractionContacts(const AtomicGroup& system,
                               const vGroup& probes,
                               const vGroup& targets,
                               const double inner_radius,
                               const double outer_radius,
                               const bool symmetry) 
{

    FContactsList fclist;
    
    for (uint j=0; j<probes.size(); ++j) {
        AtomicGroup nearby = pickNearbyAtoms(probes[j], system, outer_radius, symmetry);
        vector<double> v = fractionContactsToProbe(probes[j], nearby, targets, inner_radius, outer_radius, symmetry);
        fclist.push_back(v);
    }

    return(fclist);
}



vector<double> average(const FContactsList& f) 
{
    vector<double> avgs(f[0].size(), 0.0);
    for (uint i=0; i<avgs.size(); ++i) {
        for (uint j=0; j<f.size(); ++j)
            avgs[i] += f[j][i];
        avgs[i] /= f.size();
    }

    return(avgs);
}




int main(int argc, char *argv[]) {
    string hdr = invocationHeader(argc, argv);

    opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
    opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices();
    ToolOptions* topts = new ToolOptions;

    opts::AggregateOptions options;
    options.add(bopts).add(tropts).add(topts);

    if (!options.parse(argc, argv))
        exit(-1);

    AtomicGroup model = tropts->model;
    pTraj traj = tropts->trajectory;
    vector<uint> indices = tropts->frameList();

    AtomicGroup probe = selectAtoms(model, topts->probe_selection);
    AtomicGroup system = selectAtoms(model, topts->system_selection);
    

    // Build each of the requested targets...
    vGroup targets;
    for (vector<string>::iterator i = topts->target_selections.begin(); i != topts->target_selections.end(); ++i)
        targets.push_back(selectAtoms(model, *i));


    // If splitting, then split based on presence of connectivity...
    vGroup myselves;
    if (topts->auto_split) {
        if (probe.hasBonds())
            myselves = probe.splitByMolecule();
        else
            myselves = probe.splitByUniqueSegid();
        
    } else
        myselves.push_back(probe);


    // Size of the output matrix
    uint rows = indices.size();
    uint cols = targets.size() + 1;

    uint t = 0;
    DoubleMatrix M(rows, cols);

    // Setup our progress counter since this can be a time-consuming
    // program, but only if verbose output is requested.
    PercentProgressWithTime watcher;
    ProgressCounter<PercentTrigger, EstimatingCounter> slayer(PercentTrigger(0.1), EstimatingCounter(indices.size()));
    slayer.attach(&watcher);
    if (bopts->verbosity)
        slayer.start();

    for (vector<uint>::iterator frame = indices.begin(); frame != indices.end(); ++frame) {
        traj->readFrame(*frame);
        traj->updateGroupCoords(model);

        if (topts->symmetry && !model.isPeriodic()) {
            cerr << "ERROR - the trajectory must be periodic to use --reimage\n";
            exit(-1);
        }

        M(t, 0) = t;

        FContactsList fcl = fractionContacts(system, myselves, targets, topts->inner_cutoff, topts->outer_cutoff, topts->symmetry);
        vector<double> avg = average(fcl);
        for (uint i=0; i<avg.size(); ++i)
            M(t, i+1) = avg[i];
        
        ++t;
        if (bopts->verbosity)
            slayer.update();
    }

    if (bopts->verbosity)
        slayer.finish();

    writeAsciiMatrix(cout, M, hdr);
}
