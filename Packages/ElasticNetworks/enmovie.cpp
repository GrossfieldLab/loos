/*
  enmovie

  Elastic Network MOde VIsualizEr


  Usage:
    enmovie [options] model-name eigenvector-matrix output-name

  Notes:
    use the "--help" option for more information about how to run...

*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008 Tod D. Romo
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
#include <iterator>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include <cassert>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


typedef Math::Matrix<double> Matrix;




// --------------------------------------------------------------------------------

string vec_name;
vector<int> cols;
vector<double> scales;
double global_scale = 1.0;
bool uniform;
bool invert = false;
string model_name;
string map_name;
string selection;
bool autoscale, square;
double autolength;
string svals_file;
int verbosity = 0;
int offset = 0;
int nsteps = 100;
string supersel;
bool tag = false;



// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("mode,M", po::value< vector<string> >(&strings), "Modes to use")
      ("autoscale,A", po::value<bool>(&autoscale)->default_value(true), "Automatically scale vectors")
      ("autolength,L", po::value<double>(&autolength)->default_value(2.0), "Length of average vector in Angstroms")
      ("svals,S", po::value<string>(&svals_file), "Scale columns by singular values from file")
      ("pca", "Vectors are from PCA (sets square=1, invert=0, offset=0)")
      ("enm", "Vectors are from ENM (sets square=0, invert=1, offset=6)")
      ("superset,U", po::value<string>(&supersel)->default_value("all"), "Superset to use for frames in the output")
      ("tag,T", po::value<bool>(&tag)->default_value(false), "Tag ENM atoms with 'E' alt-loc")
      ("square", po::value<bool>(&square)->default_value(true), "square the singular values")
      ("invert", po::value<bool>(&invert)->default_value(false), "Invert singular values (ENM)")
      ("scale", po::value< vector<double> >(&scales), "Scale the requested columns")
      ("global", po::value<double>(&global_scale)->default_value(1.0), "Global scaling")
      ("uniform", po::value<bool>(&uniform)->default_value(false), "Scale all elements uniformly")
      ("map", po::value<string>(&map_name), "Use a map file to map LSV/eigenvectors to atomids")
      ("offset", po::value<int>(&offset), "Added to mode indices to select columns in eigenvector matrix");
  }

  bool postConditions(po::variables_map& vm) {

    if (vm.count("enm")) {
      square = false;
      invert = true;
      offset = 6;
    } else if (vm.count("pca")) {
      square = true;
      invert = false;
      offset = 0;
    }

    if (strings.empty())
      cols.push_back(0);
    else
      cols = parseRangeList<int>(strings);

    for (uint i=0; i<cols.size(); ++i)
      cols[i] += offset;

    if (!scales.empty()) {
      if (scales.size() != cols.size()) {
        cerr << "ERROR - You must have the same number of scalings as columns or rely on the global scaling\n";
        return(false);
      }
    } else {
      for (uint i=0; i<cols.size(); ++i)
        scales.push_back(1.0);
    }


    return(true);
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("columns='%s', global=%f, uniform=%d, map='%s', autoscale=%d, autolength=%f, svals='%s', square=%d, invert=%d, offset=%d, tag=%d, superset='%s'")
      % vectorAsStringWithCommas<string>(strings)
      % global_scale
      % uniform
      % map_name
      % autoscale
      % autolength
      % svals_file 
      % square
      % invert
      % offset
      % tag
      % supersel;
    
    oss << "scale='";
    copy(scales.begin(), scales.end(), ostream_iterator<double>(oss, ","));
    oss << "'";
    return(oss.str());
  }

  vector<string> strings;
};
// @endcond


string fullHelpMessage(void){
  string msg =                                                                                  
    "\n"                                                                                      
    "SYNOPSIS\n"                                                                      
    "\n"  
    "\tCreate a representation of motion along the mode(s) of an ENM\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "It is often informative to visualize the modes of motion predicted\n"
    "by an ENM in addition to plotting eigenvectors.  enmovie creates a dcd\n"
    "and an accompanying pdb for this purpose.  A 100 frame trajectory is \n"
    "made and the beads follow a given eigenvector(s).\n"
     "\n"
    "* PCA vs ENM *\n"
    "Enmovie should use different options depending on whether the eigenvectors come\n"
    "from a PCA or an ENM.  The --enm and --pca flags configure porcupine to expect\n"
    "the appropriate input.  If neither flag is given, then PCA is assumed.\n"
    "For PCA results, the first mode is in the first column.  LOOS\n"
    "calculates a PCA using the singular value decomposition, so the 'eigenvalues' are\n"
    "actually singular values and need to be squared.  For typical ENMs, the first 6\n"
    "eigenvectors correspond to rigid-body motion and are zero, and hence skipped.\n"
    "In addition, the magnitude of the fluctuations are the inverse of the eigenvalues.\n"
    "\n"
    "* Scaling and Autoscaling *\n"
    "There are several different ways the individual vectors can be scaled.  The default\n"
    "is to automatically determine a scaling such that the largest average displacement\n"
    "is 2 Angstroms.  If multiple modes are being used, then the corresponding eigenvector\n"
    "can be used so the relative lengths are correct.  When used with autoscaling, the\n"
    "the relative lengths are maintained.  In addition, an explicit scaling can be used\n"
    "for each mode.  If autoscaling or eigenvectors are used, then this is applied -after-\n"
    "both of those.  Finally, a global scaling can be applied.  To see the scaling used\n"
    "turn on verbose output (-v1).  For more details about exactly what scaling is used,\n"
    "set verbosity greater than 1 (-v2).\n"
    "\n"
    "In general, the default options should be fine for visualization.  If you are using\n"
    "more than one mode, then include the eigenvectors to preserve the relative scalings\n"
    "between the modes.\n"
    "\n"
    "* Supersets *\n"
    "Some visualization programs require more atoms than what the PCA/ENM used in order\n"
    "to get the structure correct (such as ribbons representations).  Including all atoms\n"
    "can solve this problem.  Alternatively, sometimes extra atoms are required to provide\n"
    "context to the region of interest, such as the extracellular loops in GPCRs.  You can\n"
    "control what atoms are written to the trajectory with the superset selection.  This\n"
    "lets you add back in atoms that were excluded by the PCA/ENM.  The catch is that they\n"
    "will not move in the trajectory, resulting in distorted bonds/connections.  The default\n"
    "is to include all atoms in the output.  If you want only the PCA/ENM region, then use\n"
    "the same selection for the superset as the vector selection.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tenmovie model.pdb pca_U.asc\n"
    "This example uses the first mode, assumes a PCA result,\n"
    "and autoscales the vectors.  Creates output.pdb and output.dcd and\n"
    "the trajectory has 100 frames.\n"
    "\n"
    "\tenmovie --pca -S pca_s.asc -M 0:3 -p modes model.pdb pca_U.asc\n"
    "This example again uses the first three modes, autoscales, and also\n"
    "scales each mode by the corresponding singular value.  It explicitly uses\n"
    "a PCA result.  It creates modes.pdb and modes.dcd with 100 frames.\n"
    "\n"
    "\tenmovie --enm -S enm_s.asc -M 0:3 -p modes model.pdb enm_U.asc\n"
    "This example is the same as above, but expects an ENM result (inverting the\n"
    "eigenvalues, and skipping the first 6 eigenpairs.\n"
    "\n"
    "\tenmovie -S pca_s.asc -M 0,3,7 -L 3 -p modes model.pdb pca_U.asc\n"
    "A PCA result is assumed, the first, fourth, and eighth mode are used, autoscaling\n"
    "is turned on with a length of 3 Angstroms.  The singular values are also included.\n"
    "The output prefix is modes."
    "\n"
    "\tenmovie --enm -S enm_s.asc -M 0,1 -A 0 -p modes --global 50 model.pdb enm_U.asc\n"
    "An ENM result is expected and the first two modes are used.  Autoscaling is disabled.\n"
    "Each mode is scaled by the corresponding eigenvalue (inverted, since this is an ENM).\n"
    "A global scaling of 50 is applied to all modes.\n"
    "\n"
    "\tenmovie --pca -S pca_s.asc -M 0:3 -p modes -U 'name == \"CA\"' model.pdb pca_U.asc\n"
    "This example again uses the first three modes, autoscales, and also\n"
    "scales each mode by the corresponding singular value.  It explicitly uses\n"
    "a PCA result.  It creates modes.pdb and modes.dcd with 100 frames.\n"
    "The default selection is to use CAs for the eigenvectors, and the -U option\n"
    "causes the output trajectory to only include CAs.\n"
    "\n"
    "\tenmovie -p pca_mode1 -t xtc model.pdb pca_U.asc\n"
    "This example uses the first mode, assumes a PCA result,\n"
    "and autoscales the vectors.  Creates pca_mode1.pdb and pca_mode1.xtc and\n"
    "the trajectory has 100 frames.\n"
    "\n"
    "SEE ALSO\n"
    "\n"
    "\tsvd, big-svd, svdcolmap, anm, gnm, vsa, porcupine\n"
    ;

  return(msg);
}


string generateSegid(const uint n) {
  stringstream s;

  s << boost::format("P%03d") % n;
  return(s.str());
}



// Accept a map file mapping the vectors (3-tuples in the rows) back
// onto the appropriate atoms

vector<int> readMap(const string& name) {
  ifstream ifs(name.c_str());
  if (!ifs) {
      cerr << "Error- cannot open " << name << endl;
      exit(-1);
  }
  
      
  string line;
  uint lineno = 0;
  vector<int> atomids;

  while (!getline(ifs, line).eof()) {
    stringstream ss(line);
    int a, b;
    if (!(ss >> a >> b)) {
      cerr << "ERROR - cannot parse map at line " << lineno << " of file " << name << endl;
      exit(-10);
    }
    atomids.push_back(b);
  }

  return(atomids);
}


// Fake the mapping, i.e. each vector corresponds to each atom...

vector<int> fakeMap(const AtomicGroup& g) {
  AtomicGroup::const_iterator ci;
  vector<int> atomids;

  for (ci = g.begin(); ci != g.end(); ++ci)
    atomids.push_back((*ci)->id());

  return(atomids);
}


// Record the atomid's for each atom in the selected subset.  This
// allows use to map vectors back onto the correct atoms when they
// were computed from a subset...

vector<int> inferMap(const AtomicGroup& g, const string& sel) {
  AtomicGroup subset = selectAtoms(g, sel);
  vector<int> atomids;

  AtomicGroup::iterator ci;
  for (ci = subset.begin(); ci != subset.end(); ++ci)
    atomids.push_back((*ci)->id());

  
  return(atomids);
}


double subvectorSize(const Matrix& U, const vector<double>& scaling, const uint j) {
  
  GCoord c;
  for (uint i=0; i<cols.size(); ++i) {
    GCoord v(U(j, i), U(j+1, i), U(j+2, i));
    c += scaling[i] * v;
  }

  return(c.length());
}


vector<double> determineScaling(const Matrix& U) {
  
  uint n = cols.size();
  vector<double> scaling(n, 1.0);
  vector<double> svals(n, 1.0);

  // First, handle singular values, if given
  if (!svals_file.empty()) {
    Matrix S;
    readAsciiMatrix(svals_file, S);
    if (verbosity > 1)
      cerr << "Read singular values from file " << svals_file << endl;
    if (S.cols() != 1) {
      cerr << boost::format("Error- singular value file is %d x %d, but it should be a %d x 1\n") % S.rows() % S.cols() % U.rows();
      exit(-2);
    }
    for (uint i = 0; i<n; ++i) {
      scaling[i] = S[cols[i]];
      if (square)
        scaling[i] *= scaling[i];
      if (invert && scaling[i] != 0.0)
        scaling[i] = 1.0 / scaling[i];
      svals[i] = scaling[i];
    }
  }

  double avg = 0.0;
  if (autoscale) {
    for (uint i=0; i<U.rows(); i += 3)
      avg += subvectorSize(U, scaling, i);
    avg /= U.rows()/3.0;

    for (uint i=0; i<n; ++i)
      scaling[i] *= autolength / avg;
  }

  // Incorporate additional scaling...
  if (verbosity > 1) {
    cerr << "Average subvector size was " << avg << endl;
    cerr << boost::format("%4s %4s %15s %15s\n") % "col" % "mode" % "sval" % "scale";
    cerr << boost::format("%4s %4s %15s %15s\n") % "----" % "----" % "---------------" % "---------------";
  }
  for (uint i=0; i<n; ++i) {
    scaling[i] *= scales[i] * global_scale;
    if (verbosity > 1)
      cerr << boost::format("%4d %4d %15.5f %15.5f\n") % cols[i] % (cols[i]-offset) % svals[i] % scaling[i];
    else if (verbosity > 0)
      cerr << boost::format("Scaling column %d by %f\n") % cols[i] % scaling[i];
  }

  return(scaling);
}


AtomicGroup renumberAndMapBonds(const AtomicGroup& model, const AtomicGroup& subset) {

  AtomicGroup renumbered = subset.copy();

  if (!renumbered.hasBonds()) {
    renumbered.renumber();
    return(renumbered);
  }

  AtomicGroup sorted = model.copy();
  sorted.sort();

  vector<int> idmap(model.size());
  for (uint i=0; i<renumbered.size(); ++i)
    idmap[renumbered[i]->index()] = i;
  
  renumbered.pruneBonds();
  for (uint i=0; i<renumbered.size(); ++i) {
    if (! renumbered[i]->hasBonds())
      continue;
    vector<int> bondlist = renumbered[i]->getBonds();
    vector<int> newbonds(bondlist.size());
    for (uint j=0; j<bondlist.size(); ++j) {
      pAtom a = sorted.findById(bondlist[j]);
      if (a == 0) {
	cerr << "Error- could not find atom id " << bondlist[j] << " in model\n";
	exit(-2);
      }
      newbonds[j] = idmap[a->index()]+1;
    }
    renumbered[i]->setBonds(newbonds); 
  }
  renumbered.renumber();
  return(renumbered);
}


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::OutputPrefix* popts = new opts::OutputPrefix("output");
  opts::OutputTrajectoryTypeOptions* otopts = new opts::OutputTrajectoryTypeOptions;
  opts::BasicSelection* sopts = new opts::BasicSelection("name == 'CA'");
  opts::ModelWithCoords* mopts = new opts::ModelWithCoords;
  ToolOptions* topts = new ToolOptions;
  opts::RequiredArguments *ropts = new opts::RequiredArguments;
  ropts->addArgument("lsv", "left-singular-vector-file");

  opts::AggregateOptions options;
  options.add(bopts).add(popts).add(otopts).add(sopts).add(mopts).add(topts).add(ropts);
  if (!options.parse(argc, argv)) {
    cerr << "***WARNING***\n";
    cerr << "The interface to enmovie has changed significantly\n";
    cerr << "and is not compatible with previous versions.  See the\n";
    cerr << "help info above, or the --fullhelp guide.\n";
    exit(-1);
  }

  verbosity = bopts->verbosity;

  // First, read in the LSVs
  Matrix U;
  readAsciiMatrix(ropts->value("lsv"), U);
  uint m = U.rows();

  vector<double> scalings = determineScaling(U);

  // Read in the average structure...
  AtomicGroup avg = mopts->model;
  AtomicGroup superset = selectAtoms(mopts->model, supersel);

  vector<int> atomids;
  if (map_name.empty()) {
    if (sopts->selection.empty())
      atomids = fakeMap(avg);
    else
      atomids = inferMap(avg, sopts->selection);

  } else
    atomids = readMap(map_name);
  
  // Double check size of atomid map
  if (atomids.size() * 3 != m) {
    cerr << boost::format("Error - The vector-to-atom map (provided or inferred) has %d atoms, but expected %d.\n") %
      atomids.size() % (m / 3);
    exit(-1);
  }

  

  pTrajectoryWriter traj = otopts->createTrajectory(popts->prefix);

  // We'll step along the Eigenvectors using a sin wave as a final scaling...
  double delta = 2*PI/nsteps;

  // Loop over requested number of frames...
  for (int frameno=0; frameno<nsteps; frameno++) {
    double k = sin(delta * frameno);

    // Have to make a copy of the atoms since we're computing a
    // displacement from the model structure...
    AtomicGroup frame = superset.copy();
    AtomicGroup frame_subset = selectAtoms(frame, sopts->selection);

    // Loop over all requested modes...
    for (uint j = 0; j<cols.size(); ++j) {
      // Loop over all atoms...
      for (uint i=0; i<frame_subset.size(); i++) {
	GCoord c = frame_subset[i]->coords();

	// This gets the displacement vector for the ith atom of the
	// ci'th mode...
        GCoord d( U(i*3, cols[j]), U(i*3+1, cols[j]), U(i*3+2, cols[j]) );
	if (uniform)
	  d /= d.length();
	
	c += k * scalings[j] * d;

	// Stuff those coords back into the Atom object...
	frame_subset[i]->coords(c);
        if (tag)
          frame_subset[i]->chainId("E");
      }
    }

    if (k == 0) {
      // Write out the selection, converting it to a PDB
      string outpdb(popts->prefix + ".pdb");
      ofstream ofs(outpdb.c_str());
      PDB pdb;
      AtomicGroup structure = renumberAndMapBonds(avg, frame);

      pdb = PDB::fromAtomicGroup(structure);
      pdb.remarks().add(hdr);
      pdb.renumber();
      ofs << pdb;
    }


    // Now that we've displaced the frame, add it to the growing trajectory...
    traj->writeFrame(frame);
  }
}
