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
        "Determine what fraction of contacts with a probe belong to the specified targets\n"
        "\n"
        "DESCRIPTION\n"
        "\tfcontacts can be used to classify the fraction of contacts made with\n"
        "a probe selection as belonging to different target selections.\n"
        "For each atom in the probe selection, all atoms that are within outer radius\n"
        "but farther than inner radius are counted as contacts.  The number of contacting\n"
        "atoms that belong to each target selection are then counted and the corresponding\n"
        "fraction is written out.  This is repeated for each time step, so the final\n"
        "output is a matrix with time increasing along the rows and each column is the\n"
        "fractional contact corresponding to the targets listed on the command line.\n"
        "\n"
        "\tThe autoself option splits the probe selection based on connectivity (if present)\n"
        "or segid.  The algorithm will iterate over each molecule and report the average\n"
        "fractional contact.\n"
        "\n"
        "\tThe exclude option determines whether the entire probe molecule is excluded\n"
        "from the total contacts or not.  For example, if the probe is a side-chain\n"
        "of a peptide and the outer radius is set long, then there will be additional\n"
        "peptide contacts that you may not want to be included in the contact count.\n"
        "Turning the exclude option on will tell fcontacts to find all atoms connected\n"
        "to the probe selection (after splitting, if so requested) and ignore these atoms\n"
        "when calculating the total number of contacts.\n"
        "\n"
        "The selection option in fcontacts is used as a 'pre-filter' for all subsequent\n"
        "selections.  This is useful for excluding hydrogens, for example.\n"
        "\n"
        "\tfcontacts --inner=0 --outer=4.5 model.psf traj.dcd 'segid == \"PEPT\"'\\\n"
        "\t          'resname == \"PEGL\"' 'segid == \"BULK\"'\n"
        "This example counts contacts as any atom with 4.5 angstroms and prints out\n"
        "the fraction of contacts with atoms having segid PEPT and PEGL residues vs\n"
        "bulk water (segid BULK)\n"
        "\n"
        "fcontacts --inner=2 --outer=5 model.psf traj.dcd 'resname == \"TRP\"'\\\n"
        "\t        'resname == \"PEGL\"' 'resname == \"PGGL\"' 'segid == \"BULK\"'\n"
        "This example counts contacts as any atom-atom distance greater than 2 angstroms\n"
        "and less than 5 angstroms.  Contacts are made between tryptophan residues\n"
        "and the fraction made to PEGL vs PGGL residues and bulk water are printed out.\n"
        "\n"
        "\tfcontacts --selection '!hydrogen' model.psf traj.dcd 'segid =~ \"PE..\"'\\\n"
        "\t          'resname == \"PCGL\"' 'segid == \"BULK\"'\n"
        "This example considers ONLY heavy atoms.  Contacts are within 4 angstroms (the defaults)\n"
        "and for any atom with a segid matching the pattern PExx (e.g. PE00, PE01, PE02, ...).\n"
        "The fraction of contacts made with PCGL residues and bulk solvent are printed\n"
        ;
    
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
        auto_split(true),
        exclude_self(true),
        report_stddev(false)
        { }


    void addGeneric(po::options_description& o) {
        o.add_options()
            ("inner", po::value<double>(&inner_cutoff)->default_value(inner_cutoff), "Inner cutoff (ignore atoms closer than this)")
            ("outer", po::value<double>(&outer_cutoff)->default_value(outer_cutoff), "Outer cutoff (ignore atoms further away than this)")
            ("reimage", po::value<bool>(&symmetry)->default_value(symmetry), "Consider symmetry when computing distances")
            ("split", po::value<bool>(&auto_split)->default_value(auto_split), "Automatically split probe selection")
            ("exclude", po::value<bool>(&exclude_self)->default_value(exclude_self), "Exclude self from contacts")
            ("pad", po::value<double>(&pad)->default_value(pad), "Padding for filtering nearby atoms")
            ("stddev", po::value<bool>(&report_stddev)->default_value(report_stddev), "Include stddev in output");
    }

    void addHidden(po::options_description& o) {
        o.add_options()
            ("probe", po::value<string>(&probe_selection), "Probe selection")
            ("target", po::value< vector<string> >(&target_selections), "Target selections");
    }

    void addPositional(po::positional_options_description& p) {
        p.add("probe", 1);
        p.add("target", -1);
    }

    bool check(po::variables_map& map) {
        if (target_selections.empty() || probe_selection.empty())
            return(true);

        return(false);
    }


    string help() const {
        return("probe target [target ...]");
    }

    string print() const {
        ostringstream oss;
        oss << boost::format("inner=%f,outer=%f,reimage=%d,autosplit=%d,pad=%f,stddev=%d,probe='%s',targets=")
            % inner_cutoff
            % outer_cutoff
            % symmetry
            % auto_split
            % pad
            % report_stddev
            % probe_selection;

        for (uint i=0; i<target_selections.size(); ++i)
            oss << "'" << target_selections[i] << "'" << (i == target_selections.size() - 1 ? "" : ",");

        return(oss.str());
    }


    double inner_cutoff, outer_cutoff, pad;
    string probe_selection;
    bool symmetry, auto_split, exclude_self, report_stddev;
    vector<string> target_selections;
};



// Check for atom equality only through atomid...

struct IdEquals
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

    vector<double> fracts(targets.size(), 0.0);
    for (uint i=0; i<targets.size(); ++i) {
        AtomicGroup target_nearby = nearby_contacts.intersect(targets[i], IdEquals());
        if (nearby_contacts.empty())
            fracts[i] = 0.0;
        else 
            fracts[i] = static_cast<double>(target_nearby.size()) / nearby_contacts.size();
    }
    return(fracts);
}

        
FContactsList fractionContacts(const AtomicGroup& system,
                               const vGroup& probes,
                               const vGroup& excludeds,
                               const vGroup& targets,
                               const double inner_radius,
                               const double outer_radius,
                               const bool symmetry) 
{

    FContactsList fclist;
    
    for (uint j=0; j<probes.size(); ++j) {
        AtomicGroup nearby = pickNearbyAtoms(probes[j], excludeds[j], outer_radius, symmetry);
        vector<double> v = fractionContactsToProbe(probes[j], nearby, targets,
                                                   inner_radius, outer_radius, symmetry);
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


vector<double> stddevs(const FContactsList& f, const vector<double>& avgs) 
{
    vector<double> stds(f[0].size(), 0.0);
    for (uint i=0; i<stds.size(); ++i) {
        for (uint j=0; j<f.size(); ++j) {
            double d = f[j][i] - avgs[i];
            stds[i] += d*d;
        }
        stds[i] /= (f.size()-1);
    }

    return(stds);
}

        

int main(int argc, char *argv[]) {
    string hdr = invocationHeader(argc, argv);

    opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
    opts::BasicSelection* sopts = new opts::BasicSelection();
    opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices();
    ToolOptions* topts = new ToolOptions;

    opts::AggregateOptions options;
    options.add(bopts).add(sopts).add(tropts).add(topts);

    if (!options.parse(argc, argv))
        exit(-1);

    AtomicGroup model = tropts->model;
    pTraj traj = tropts->trajectory;
    vector<uint> indices = tropts->frameList();

    AtomicGroup system = selectAtoms(model, sopts->selection);
    AtomicGroup probe = selectAtoms(system, topts->probe_selection);
    

    // Build each of the requested targets...
    vGroup targets;
    for (vector<string>::iterator i = topts->target_selections.begin(); i != topts->target_selections.end(); ++i)
        targets.push_back(selectAtoms(system, *i));


    // If splitting, then split based on presence of connectivity...
    vGroup myselves;
    vGroup excludes;
    
    if (topts->auto_split) {
        if (probe.hasBonds())
            myselves = probe.splitByMolecule();
        else
            myselves = probe.splitByUniqueSegid();
        
    } else
        myselves.push_back(probe);


    if (topts->exclude_self) {
        vGroup molecules;
        if (system.hasBonds())
            molecules = system.splitByMolecule();
        else
            molecules = system.splitByUniqueSegid();

        for (vGroup::iterator i = myselves.begin(); i != myselves.end(); ++i) {
            AtomicGroup exclusive;
            for (vGroup::iterator j = molecules.begin(); j != molecules.end(); ++j)
                if (i->containsAny(*j, IdEquals()))
                    exclusive.append(*j);
            excludes.push_back(exclusive);
        }
    } else
        excludes = myselves;

    // This is system excluding requested probe atoms...
    vGroup excludeds;
    for (vGroup::iterator i = excludes.begin(); i != excludes.end(); ++i) {
        AtomicGroup pruned = system;
        pruned.remove(*i);
        excludeds.push_back(pruned);
    }
    
    
    // Size of the output matrix
    uint rows = indices.size();
    uint cols = targets.size();
    if (topts->report_stddev)
        cols *= 2;
    ++cols;
    

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

        M(t, 0) = *frame;

        FContactsList fcl = fractionContacts(system, myselves, excludeds, targets, topts->inner_cutoff, topts->outer_cutoff, topts->symmetry);
        vector<double> avg = average(fcl);
        if (topts->report_stddev) {
            vector<double> stds = stddevs(fcl, avg);
            uint k = 0;
            for (uint i=0; i<avg.size(); ++i) {
                M(t, ++k) = avg[i];
                M(t, ++k) = stds[i];
            }

        } else {

            for (uint i=0; i<avg.size(); ++i)
                M(t, i+1) = avg[i];

        }
        
        ++t;
        if (bopts->verbosity)
            slayer.update();
    }

    if (bopts->verbosity)
        slayer.finish();

    writeAsciiMatrix(cout, M, hdr);
}
