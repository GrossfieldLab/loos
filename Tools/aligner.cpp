/*
  aligner.cpp


  Aligns structures in a trajectory...

  Usage:

    aligner [options] structure-file trajectory-file output-prefix

  Notes:

  Takes two selections.  The first is the subset of atoms that will
  be used for the alignment.  The second is the subset of atoms that
  will then be transformed by the alignment and written out.  The
  alignment scheme is to calculate an average structure, align all
  frames against this average, then compute a new average.  This is
  iterated until the difference in average structures is below the
  specified tolerance.

  The output will ALWAYS be a DCD!

  Aligner will cache the entire alignment selection in memory, so
  beware potential memory issues...

*/




/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
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


string model_name, traj_name, prefix;
string alignment_string;
string transform_string;

string reference_name;
string reference_sel;

greal alignment_tol;

int maxiter;
bool center;



// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() : alignment_string("name == 'CA'"),
                  transform_string("all"),
                  reference_name(""),
                  reference_sel(""),
                  alignment_tol(1e-6),
                  maxiter(5000),
                  center(true) { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("align", po::value<string>(&alignment_string)->default_value(alignment_string), "Align using this selection")
      ("transform", po::value<string>(&transform_string)->default_value(transform_string), "Transform using this selection")
      ("maxiter", po::value<uint>(&maxiter)->default_value(maxiter), "Maximum number of iterations for alignment algorith")
      ("tolerance", po::value<double>(&alignment_tol)->default_value(alignment_tol), "Tolerance for alignment convergence")
      ("center", po::value<bool>(&center)->default_value(center), "Auto-center the trajectory using the alignment subset")
      ("reference", po::value<string>(&reference_name), "Align to a reference structure (non-iterative")
      ("refsel", po::value<string>(&reference_sel), "Selection to align against in reference (default is same as --align)");
  }

  string alignment_string, transform_string;
  string reference_name, reference_sel;
  double alignment_tol;
  uint maxiter;
  bool center;
};


// @endcond



void centerFrame(AtomicGroup& src, AtomicGroup& trg) {
  GCoord c = src.centroid();
  trg.translate(-c);
}

void savePDB(const string& fname, const string& meta, const AtomicGroup& grp) {
  AtomicGroup dup = grp.copy();
  PDB pdb = PDB::fromAtomicGroup(dup);
  pdb.pruneBonds();
  pdb.renumber();
  pdb.remarks().add(meta);
  ofstream ofs(fname.c_str());
  ofs << pdb;
}




int main(int argc, char *argv[]) {

  // Parse command-line options, cache invocation header for later use...
  string header = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::OutputPrefix* prefopts = new opts::OutputPrefix;
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  ToolOptions* toolopts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(prefopts).add(tropts).add(toolopts);
  if (!options.parse(argc, argv))
    exit(-1);

  // Read the inputs...
  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;

  // Get the selections (subsets) to operate over
  AtomicGroup align_sub = selectAtoms(model, toolopts->alignment_string);

  AtomicGroup applyto_sub = selectAtoms(model, toolopts->transform_string);

  // Now do the alignin'...
  unsigned int nframes = traj->nframes();

  
  if (toolopts->reference_name.empty()) {

    // Read in the trajectory frames and extract the coordinates for the
    // aligning subset...
    vector<AtomicGroup> frames;
    while (traj->readFrame()) {
      traj->updateGroupCoords(align_sub);
      AtomicGroup subcopy = align_sub.copy();
      frames.push_back(subcopy);
    }
    
    boost::tuple<vector<XForm>,greal, int> res = iterativeAlignment(frames, toolopts->alignment_tol, toolopts->maxiter);
    greal final_rmsd = boost::get<1>(res);
    cerr << "Final RMSD between average structures is " << final_rmsd << endl;
    cerr << "Total iters = " << boost::get<2>(res) << endl;
    
    vector<XForm> xforms = boost::get<0>(res);
    
    // Zzzzap our stored groups...
    frames.clear();
    
    // Setup for writing DCD...
    DCDWriter dcdout(prefopts->prefix + ".dcd");
    dcdout.setHeader(applyto_sub.size(), nframes, 1e-3, traj->hasPeriodicBox());
    dcdout.setTitle(header);
    dcdout.writeHeader();
    
    // Now apply the alignment transformations to the requested subsets
    for (unsigned int i = 0; i<nframes; i++) {
      traj->readFrame(i);
      traj->updateGroupCoords(model);
      applyto_sub.applyTransform(xforms[i]);
      
      if (center)
        centerFrame(align_sub, applyto_sub);
      dcdout.writeFrame(applyto_sub);
      
      if (i == 0) 
        savePDB(prefopts->prefix + ".pdb", header, applyto_sub);
    }
    
  } else {
    
    AtomicGroup reference = createSystem(toolopts->reference_name);
    
    string refsel = toolopts->reference_sel.empty() ? toolopts->alignment_string : toolopts->reference_sel;
    AtomicGroup refsub = selectAtoms(reference, refsel);

    if (refsub.size() != align_sub.size()) {
      cerr << boost::format("ERROR- alignment subset has %d atoms but reference subset has %d.  They must match.\n") % align_sub.size() % refsub.size();
      exit(-10);
    }

    DCDWriter dcdout(prefopts->prefix + ".dcd");
    dcdout.setHeader(applyto_sub.size(), nframes, 1e-3, traj->hasPeriodicBox());
    dcdout.setTitle(header);
    dcdout.writeHeader();


    bool first = true;
    while (traj->readFrame()) {
      traj->updateGroupCoords(model);
      GMatrix M = align_sub.superposition(refsub);
      XForm W(M);
      applyto_sub.applyTransform(W);

      if (center)
        centerFrame(align_sub, applyto_sub);
      dcdout.writeFrame(applyto_sub);

      if (first) {
        savePDB(prefopts->prefix + ".pdb", header, applyto_sub);
        first = false;
      }
    }

  }
}



    
