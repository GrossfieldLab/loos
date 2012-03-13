/*
  water-inside.cpp


  (c) 2008 Tod D. Romo,
      Grossfield Lab,
      University of Rochester Medical and Dental School


  Applies a given set of criteria to determine whether or not a water
  is inside a protein.  A matrix is then built up where each column
  represents a timepoint in the trajectory and each row is the
  internal water state (i.e. 1 = water is inside, 0 = water is not
  inside)

  Also tracks the volume of the probe region (i.e. what's defined as
  inside, if possible) and writes out a list of atomids that describe
  which atoms go with which rows of the matrix.

*/

#include <loos.hpp>



#include <loos.hpp>
#include <cmath>
#include <boost/format.hpp>

#include <DensityGrid.hpp>
#include <DensityTools.hpp>
#include <DensityOptions.hpp>

using namespace std;
using namespace loos;
using namespace loos::DensityTools;



typedef Math::Matrix<int, Math::ColMajor>    Matrix;


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "\tClassify waters as inside a protein or not in a trajectory\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\twater-inside applies a user-specified set of criteria to determine\n"
    "whether or not water is inside a protein.  A matrix is constructed where\n"
    "each column is a time-series for each water.  A 1 means the corresponding\n"
    "water is inside the protein, and 0 means it's not. The volume of the probe region\n"
    "is also tracked (if possible), and written out separately.  In addition, a file\n"
    "is written that maps the water atomids to the columns of the matrix.\n"
    "See water-hist for more information about internal-water criteria.\n"
    "\n"
    "\nEXAMPLES\n"
    "\twater-inside --prefix water foo.pdb foo.dcd\n"
    "This example will use the axis filter for water (i.e. water atoms within\n"
    "the default radius of 10 Angstroms from the first principal axis of the protein\n"
    "selection.  The default water selection (name == 'OH2') and protein selection\n"
    "(name == 'CA') are used.  The output prefix is set to 'water', so 'water.asc',\n"
    "'water.vol', and 'water.atoms' will be created containing the time-series matrix,\n"
    "the internal water region volume, and the atom mapping respectively.\n\n"
    "\twater-inside --mode radius --radius 5 --prot 'resid == 65' --prefix pocket foo.pdb foo.dcd\n"
    "This example will find water atoms (using the default selection) that are within\n"
    "5 Angstroms of any atom in residue 65, and use the output prefix 'pocket'.\n"
    "\n"
    "NOTES\n"
    "\tLOOS does not care what is called a protein or water.  You can use any selection,\n"
    "for example, to track ligands, or lipids, etc.\n"
    "\n"
    "SEE ALSO\n"
    "\twater-hist\n"
    ;

  return(msg);
}


    

void writeAtomIds(const string& fname, const AtomicGroup& grp, const string& hdr) {
  ofstream ofs(fname.c_str());
  AtomicGroup::const_iterator ci;

  uint t = 0;
  ofs << "# " << hdr << endl;
  ofs << "# i\tatomid(i)\tresidue(i)\n";
  for (ci = grp.begin(); ci != grp.end(); ++ci) {
    ofs << boost::format("%u\t%d\t%s-%d\n") % t++ % (*ci)->id() % (*ci)->name() % (*ci)->resid();
  }
}


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* basopts = new opts::BasicOptions(fullHelpMessage());
  opts::OutputPrefix* prefopts = new opts::OutputPrefix;
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  opts::BasicWater* watopts = new opts::BasicWater;

  opts::AggregateOptions options;
  options.add(basopts).add(prefopts).add(tropts).add(watopts);
  if (!options.parse(argc, argv))
    exit(-1);



  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;
  vector<uint> frames = tropts->frameList();

  AtomicGroup subset = selectAtoms(model, watopts->prot_string);
  AtomicGroup waters = selectAtoms(model, watopts->water_string);

  uint m = waters.size();
  uint n = traj->nframes();
  Matrix M(m,n);
  Math::Matrix<double> V(n, 1);
  cerr << boost::format("Water matrix is %d x %d.\n") % m % n;

  uint i = 0;
  cerr << "Processing - ";

  for (vector<uint>::iterator t = frames.begin(); t != frames.end(); ++t) {
    if (i % 100 == 0)
      cerr << ".";

    traj->readFrame(*t);
    traj->updateGroupCoords(model);

    vector<int> mask = watopts->filter_func->filter(waters, subset);
    if (mask.size() != m) {
      cerr << boost::format("ERROR - returned mask has size %u but expected %u.\n") % mask.size() % m;
      exit(-10);
    }

    for (uint j=0; j<m; ++j)
      M(j,i) = mask[j];

    V(i,0) = watopts->filter_func->volume();
    ++i;
  }

  cerr << " done\n";
  writeAsciiMatrix(prefopts->prefix + ".asc", M, hdr);
  writeAsciiMatrix(prefopts->prefix + ".vol", V, hdr);
  writeAtomIds(prefopts->prefix + ".atoms", waters, hdr);
}
