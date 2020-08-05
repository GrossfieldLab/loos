/*
  svd.cpp

  Computes the SVD for a trajectory.  Writes out the SVD as an
  OCTAVE-formatted text file.
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
#include <boost/filesystem.hpp>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

#if defined(__linux__)
extern "C" {
  void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
  void sgesvd_(char*, char*, int*, int*, float*, int*, float*, float*, int*, float*, int*, float*, int*, int*);
}
#endif


typedef double svdreal;
#define SVDFUNC  dgesvd_


typedef Math::Matrix<svdreal, Math::ColMajor> Matrix;



// Globals
string header("NO HEADER SPECIFIED");
string prefix("output");


// @cond TOOLS_INTERNAL

class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() :
    alignment_string("name == 'CA'"),
    svd_string("name == 'CA'"),
    noalign(false),
    include_source(false),
    alignment_tol(1e-6),
    splitv(true),
    autoname(true),
    terms(0)
  { }


  void addGeneric(po::options_description& o) {
    o.add_options()
      ("align,A", po::value<string>(&alignment_string)->default_value(alignment_string), "Selection to align with")
      ("svd,S", po::value<string>(&svd_string)->default_value(svd_string), "Selection to calculate the SVD of")
      ("tolerance", po::value<double>(&alignment_tol)->default_value(alignment_tol), "Tolerance for iterative alignment")
      ("noalign,N", po::value<bool>(&noalign)->default_value(noalign), "Do NOT align the frames of the trajectory")
      ("source", po::value<bool>(&include_source)->default_value(include_source), "Write out source conformation matrix")
      ("splitv", po::value<bool>(&splitv)->default_value(splitv), "Automatically split V matrix (when using multiple trajectories)")
      ("autoname", po::value<bool>(&autoname)->default_value(autoname), "Automatically name V files based on traj filename")
      ("terms", po::value<uint>(&terms), "# of terms of the SVD to output");
  }


  bool postConditions(po::variables_map& vm) {
    if (autoname)
      splitv = true;

    return(true);
  }


  
  string print() const {
    ostringstream oss;

    oss << boost::format("align='%s', svd='%s', tolerance=%f, noalign=%d, source=%d, splitv=%d, autoname=%d, terms=%d")
      % alignment_string
      % svd_string
      % noalign
      % include_source
      % alignment_tol
      % splitv
      % autoname
      % terms;
    return(oss.str());
  }


  string alignment_string, svd_string;
  bool noalign, include_source;
  double alignment_tol;
  bool splitv, autoname;
  uint terms;
};

// @endcond



string fullHelpMessage(void)
{
string s =
  "\n"
  "SYNOPSIS\n"
  "\n"
  "Calculate the principal components of a simulation using\n"
  "the singular value decomposition\n"
  "\n"
  "DESCRIPTION\n"
  "\n"
  "This tool performs a principal component analysis (PCA)\n"
  "on the trajectory.  This technique computes a new coordinate\n"
  "system such that the largest concerted motions are on the 1st\n"
  "axis (the 1st principal component).  This effectively reduces\n"
  "relevant dimensionality of the system by resolving the most\n"
  "collective motions (those with the largest covariance) followed\n"
  "by those with the 2nd largest covariance, etc...\n"
  "\n"
  "This technique is also referred to in the literature as essential\n"
  "dynamics.  This tool performs the PCA using a technique called the\n"
  "singular value decomposition (SVD).  There are several output files\n"
  "that can be used for  numerous analyses.  A list of the files and\n"
  "their contents follows.  For these descriptions assume an SVD is\n"
  "\t\t A = UsV*\n"
  "where the matrix A contains the coordinates of the atoms for every\n"
  "frame in the trajectory:\n"
  "\toutput_s.asc   - singular values (square roots of eigenvalues)\n"
  "\toutput_U.asc   - left singular vectors (lsv, direction of each PC)\n"
  "\toutput_V.asc   - right singular vectors (rsv, motion of a frame \n"
  "\t                    projected onto the PC with the same index\n" 
  "\toutput.map     - mapping of selection onto rows of output matrices\n"
  "\toutput_avg.pdb - average structure across the trajectory\n"
  "\n"
  "\n"
  "UNITS AND PCA COMPARISON\n"
  "\n"
  "The left and right singular vectors are column vectors, meaning that\n"
  "each column of the matrix (U and V respectively) is a vector.  These\n"
  "vectors must have length 1, so their elements are normalized.  The\n"
  "left singular vectors (LSVs) are normalized by 1/sqrt(L) there L is the\n"
  "length of the trajectory (i.e. number of frames).  The right singular vectors\n"
  "(RSVs) are normalized by 1/sqrt(3n) where n is the number of atoms.\n"
  "The LSVs are the same as the eigenvectors from a traditional PCA.\n"
  "The singular values are just the square roots of the PCA eigenvalues,\n"
  "and are in Angstroms.\n"
  "\n"
  //
  "EXAMPLES\n"
  "\n"
  "svd -A 'name==\"CA\"' -S 'name==\"CA\"' model.pdb traj.dcd\n"
  "\tComputes the PCA of the CA's in model.pdb across the entire trajectory\n"
  "\ttraj.dcd.  The file output_U.asc contains the LSVs, which point in the\n"
  "\tdirection of motion associated with each eigenvalue.  The square roots\n"
  "\tof the eigenvalues are contained in output_s.asc.  The \"-A\" option says\n"
  "\tthat the trajectory will be aligned using the CA's prior to the PCA. \n"
  "\tSee \"aligner\" for more details on trajectory alignment.\n"
  "\n"
  "svd -k25 -A 'name==\"CA\"' -S 'name==\"CA\"' model.pdb traj.dcd\n"
  "\tSame as the example above but here we are skipping the 1st 25 frames\n"
  "\tof the trajectory.  A common reason for this might be allowing the \n"
  "\tsystem additional sampling before data analysis.\n"
  "\n"
  "svd -r 25:5:250 -A 'name==\"CA\"' -S 'name==\"CA\"' model.pdb traj.dcd\n"
  "\tThis example uses the octave-style range info to decide which frames of\n"
  "\tthe simulation to use for the PCA. Similar to the case above we skip the\n"
  "\t1st 25 frames.  We will calculate the PCA upto frame 250, while using\n"
  "\tonly every 5th frame.  A common use for this option might be analyzing\n"
  "\tonly a specific, large feature of the simulation.\n"
  "\n"
  "svd -p svd_model -N1 -S 'segid==\"PROT\" && !(hydrogen)' model.pdb traj.dcd\n"
  "\tPerform the svd of the same simulation with a few changes.  First, the\n"
  "\toutput files have the prefix \"svd_model\" (i.e. svd_model_u.asc).  Next,\n"
  "\twe are not aligning the trajectory.  Finally, we are now computing the \n"
  "\tPCs of all heavy atoms in the protein (segid PROT).\n"
  "\t\n"
  "\t\n"
  //
  "SEE ALSO\n"
  "\n"
  "A number of LOOS analysis tools work on PCA results. Here is a partial list:\n"
  "Other programs in Tools:\n"
  "\tporcupine - Vizualization tool: create a pdb with \"sticks\" that points\n"
  "\t               in the direction of the PCs\n"
  "\tphase-pdb - Vizualization tool: use 3 PCs for vizualization\n" 
  "\tcoverlap  - calculate the covariance overlap between 2 PCA results\n"
  "\tbig-svd   - calculates svd of a trajectory too big for svd (this tool)\n"
  "\t               to handle.  Uses a slightly different algorithm\n"
  "\tkurskew   - calculates the skew and kurtosis of each column in a matrix\n" 
  "\t\n" 
  "Use with Packages/ElasticNetworks:\n"
  "\tenmovie - create a dcd of motion along a PC (for visualization)\n"
  "\t\n"
  "Also, some convergence tools make use of PCA:\n"
  "In these tools PCA is performed within the program's execution:\n"
  "\tbcom\n"
  "\tboot-bcom\n"
  "\n"
  "These tools require input of file(s) from svd:\n"
  "\trsv-coscon\n"
  "\t\n"
  "\n";

    return (s);
    }






vector<XForm> doAlign(const AtomicGroup& subset, pTraj traj, const vector<uint>& indices, const double tol) {

  boost::tuple<vector<XForm>, greal, int> res = iterativeAlignment(subset, traj, indices, tol, 100);
  vector<XForm> xforms = boost::get<0>(res);
  greal rmsd = boost::get<1>(res);
  int iters = boost::get<2>(res);

  cerr << "Subset alignment with " << subset.size()
       << " atoms converged to " << rmsd << " rmsd after "
       << iters << " iterations.\n";

  return(xforms);
}


void writeAverage(const AtomicGroup& avg) {
  PDB avgpdb = PDB::fromAtomicGroup(avg);
  avgpdb.pruneBonds();
  avgpdb.remarks().add(header);

  string fname(prefix + "_avg.pdb");
  ofstream ofs(fname.c_str());
  if (!ofs)
    throw(runtime_error("Cannot open " + fname));
  ofs << avgpdb;
}


// Calculates the transformed avg structure, then extracts the
// transformed coords from the DCD with the avg subtraced out...

Matrix extractCoords(const AtomicGroup& subset, const vector<XForm>& xforms, pTraj traj, const vector<uint>& indices) {

  AtomicGroup avg = averageStructure(subset, xforms, traj, indices);
  writeAverage(avg);

  uint natoms = subset.size();
  AtomicGroup frame = subset.copy();
  uint n = indices.size();
  uint m = natoms * 3;
  Matrix M(m, n);

  for (uint i=0; i<n; ++i) {
    traj->readFrame(indices[i]);
    traj->updateGroupCoords(frame);
    frame.applyTransform(xforms[i]);

    for (uint j=0; j<natoms; j++) {
      GCoord c = frame[j]->coords() - avg[j]->coords();
      M(j*3,i) = c.x();
      M(j*3+1,i) = c.y();
      M(j*3+2,i) = c.z();
    }
  }

  return(M);
}



void write_map(const string& fname, const AtomicGroup& grp) {
  ofstream fout(fname.c_str());

  if (!fout) {
    cerr << "Unable to open " << fname << " for output.\n";
    exit(-10);
  }

  pAtom pa;
  AtomicGroup::Iterator iter(grp);
  int i = 0;
  while ( (pa = iter()) )
    fout << i++ << "\t" << pa->id() << "\t" << pa->resid() << endl;

}


void writeMatrixChunk(opts::OutputPrefix* popts, opts::MultiTrajOptions* tropts, ToolOptions* topts, const Matrix& Vt, const Math::Range& start, const Math::Range& end, const string& header, const uint index) {
  string filename;

  if (topts->autoname) {
    boost::filesystem::path p(tropts->mtraj[index]->filename());
#if BOOST_FILESYSTEM_VERSION >= 3
    filename = p.stem().string() + "_V.asc";
#else
    filename = p.stem() + "_V.asc";
#endif
  } else {
    ostringstream oss;
    oss << boost::format("%s_V_%04d.asc") % popts->prefix % index;
    filename = oss.str();
  }

  writeAsciiMatrix(filename, Vt, header, start, end, true);
}




int main(int argc, char *argv[]) {
  header = invocationHeader(argc, argv);
  opts::BasicOptions* bhopts = new opts::BasicOptions(fullHelpMessage());
  opts::OutputPrefix* popts = new opts::OutputPrefix;
  opts::MultiTrajOptions* tropts = new opts::MultiTrajOptions;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bhopts).add(popts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  if (bhopts->verbosity)
    cout << tropts->trajectoryTable() << endl;
  
  prefix = popts->prefix;
  AtomicGroup model = tropts->model;
  pTraj ptraj = tropts->trajectory;
  vector<uint> indices = tropts->frameList();

  AtomicGroup svdsub = selectAtoms(model, topts->svd_string);

  write_map(prefix + ".map", svdsub);

  vector<XForm> xforms;
  if (topts->noalign) {
    // Make noop xforms to prevent doing any alignment...
    cerr << argv[0] << ": SKIPPING ALIGNMENT\n";
    for (uint i=0; i<indices.size(); ++i)
      xforms.push_back(XForm());
  } else {
    AtomicGroup alignsub = selectAtoms(model, topts->alignment_string);
    cerr << argv[0] << ": Aligning...\n";
    xforms = doAlign(alignsub, ptraj, indices, topts->alignment_tol);   // Honors indices
  }

  cerr << argv[0] << ": Extracting coordinates...\n";
  Matrix A = extractCoords(svdsub, xforms, ptraj, indices);   // Honors indices
  f77int m = A.rows();
  f77int n = A.cols();
  f77int sn = m<n ? m : n;


  if (topts->include_source)
    writeAsciiMatrix(prefix + "_A.asc", A, header);

  double estimate = static_cast<double>(m)*m*sizeof(svdreal) + static_cast<double>(n)*n*sizeof(svdreal) + static_cast<double>(m)*n*sizeof(svdreal) + sn*sizeof(svdreal);
  cerr << boost::format("%s: Allocating estimated %.3f GB for %d x %d SVD\n")
    % argv[0]
    % (estimate / gigabytes)
    % m
    % n;

  char jobu = 'A', jobvt = 'A';
  f77int lda = m, ldu = m, ldvt = n, lwork= -1, info;
  svdreal prework[10], *work;

  Matrix U(m,m);
  Matrix S(sn,1);
  Matrix Vt(n,n);
  
  // First, request the optimal size of the work array...
  SVDFUNC(&jobu, &jobvt, &m, &n, A.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, prework, &lwork, &info);
  if (info != 0) {
    cerr << "Error code from size request to dgesvd was " << info << endl;
    exit(-2);
  }

  lwork = (f77int)prework[0];
  estimate += lwork * sizeof(svdreal);
  cerr << argv[0] << ": SVD requests " << lwork << " extra space for a grand total of " << estimate / gigabytes << " GB\n";
  work = new svdreal[lwork];

  cerr << argv[0] << ": Calculating SVD...\n";
  Timer<WallTimer> timer;
  timer.start();
  SVDFUNC(&jobu, &jobvt, &m, &n, A.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, work, &lwork, &info);
  timer.stop();
  cerr << argv[0] << ": Done!  Calculation took " << timeAsString(timer.elapsed()) << endl;

  if (info > 0) {
    cerr << "Convergence error in dgesvd\n";
    exit(-3);
  } else if (info < 0) {
    cerr << "Error in " << info << "th argument to dgesvd\n";
    exit(-4);
  }


  Math::Range orig(0,0);
  Math::Range Usize(m,m);
  Math::Range Ssize(sn,1);
  Math::Range Vsize(sn,n);

  if (topts->terms) {
    int terms = static_cast<int>(topts->terms);
    if (terms > m || terms > sn || terms > n) {
      cerr << "ERROR- The number of terms requested exceeds matrix dimensions.\n";
      exit(-1);
    }
    Usize = Math::Range(m, terms);
    Ssize = Math::Range(terms, 1);
    Vsize = Math::Range(terms, n);
  }

  cerr << argv[0] << ": Writing results...\n";
  writeAsciiMatrix(prefix + "_U.asc", U, header, orig, Usize);
  writeAsciiMatrix(prefix + "_s.asc", S, header, orig, Ssize);

  if (topts->splitv && tropts->mtraj.size() > 1) {
    // Need to reconstruct what row-ranges correspond to the input trajectories...
    uint a = 0;
    uint curtraj = 0;
    int terms = topts->terms ? static_cast<int>(topts->terms) : sn;

    for (uint i=0; i<n; ++i) {
      MultiTrajectory::Location loc = tropts->mtraj.frameIndexToLocation(indices[i]);
      if (loc.first != curtraj) {
        writeMatrixChunk(popts, tropts, topts, Vt, Math::Range(0, a), Math::Range(terms, i), header, curtraj);
        a = i;
        curtraj = loc.first;
      }
    }

    writeMatrixChunk(popts, tropts, topts, Vt, Math::Range(0, a), Math::Range(terms, n), header, curtraj);
    
  } else
    writeAsciiMatrix(prefix + "_V.asc", Vt, header, orig, Vsize, true);
  
  cerr << argv[0] << ": done!\n";

  delete[] work;
}
