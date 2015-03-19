/*
  porcupine
  

  Creates "porcupine" plots by placing atoms at the endpoints of the
  vectors and adding a bond between them.  This is all written out a
  PDB file with CONECT records.

  Notes:

    You can use a "map" file to map the eigenvectors back onto specific
  atoms via atomid.  The eigenvector matrix is arranged so that the
  eigenvectors are stored as column-vectors, but each triplet of rows
  is the eigenvector corresponding to an individual atom.  The map
  then consists of two columns: the left is the tripled index
  (i.e. row/3) and the right column is the atomid of the corresponding
  atom.
  
     0     30
     1     42
     2     57
     3     66

  Alternatively, porcupine can infer the mapping of vectors to atoms
  by giving it the same selection you used to compute the vectors
  along with the same model...

*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo
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
#include <limits>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;



typedef Math::Matrix<float, Math::ColMajor>   Matrix;


string vec_name;
vector<int> cols;
vector<double> scales;
double global_scale = 1.0;
bool uniform;
bool invert = false;
bool double_sided;
string model_name;
string map_name;
double tip_size;
string selection;
bool autoscale, square;
double autolength;
string svals_file;
int verbosity = 0;
int offset = 0;

const string porcupine_tag("POR");
const string tip_tag("POT");

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
      ("tips,T", po::value<double>(&tip_size)->default_value(0.0), "Length (in Angstroms) to make the tip (for single-sided only)")
      ("double_sided", po::value<bool>(&double_sided)->default_value(false), "Use double-sided vectors")
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
    oss << boost::format("columns='%s', global=%f, uniform=%d, map='%s', tips=%f, double_sided=%d, autoscale=%d, autolength=%f, svals='%s', square=%d, invert=%d, offset=%d")
      % vectorAsStringWithCommas<string>(strings)
      % global_scale
      % uniform
      % map_name
      % tip_size
      % double_sided
      % autoscale
      % autolength
      % svals_file 
      % square
      % invert
      % offset;
    
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
    "\tCreate a matchstick representation of eigenvectors/left singular vectors (LSV)\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis program takes a model and a vector-matix and creates a pdb illustrating\n"
    "the direction of those vectors starting from the model structure.  \n"
    "The typical use is for illustrating the direction of motion calculated from\n"
    "a trajectory PCA or predicted from NMA of a network model.\n"
    "\n"
    "* PCA vs ENM *\n"
    "Porcupine should use different options depending on whether the eigenvectors come\n"
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
    "is to automatically determine a scaling such that the largest average drawn vector\n"
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
    "* The Model *\n"
    "The resulting PDB has the following properties...  Each mode has its own segid\n"
    "in the form 'Pnnn' there nnn is a zero-padded mode number.  Each vector has\n"
    "an atom name of 'POR' and residue name of 'POR'.  The vectors have increasing\n"
    "resids that reset for each mode.  If tips are used, then the tip atoms will\n"
    "have an atom name of 'POT'.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tporcupine model.pdb pca_U.asc >porcupine.pdb\n"
    "This example uses the first mode, assumes a PCA result,\n"
    "and autoscales the vectors.\n"
    "\n"
    "\tporcupine --pca -S pca_s.asc -M 0:3 model.pdb pca_U.asc >porcupine.pdb\n"
    "This example again uses the first three modes, autoscales, and also\n"
    "scales each mode by the corresponding singular value.  It explicitly uses\n"
    "a PCA result.\n"
    "\n"
    "\tporcupine --enm -S enm_s.asc -M 0:3 model.pdb enm_U.asc >porcupine.pdb\n"
    "This example is the same as above, but expects an ENM result (inverting the\n"
    "eigenvalues, and skipping the first 6 eigenpairs.\n"
    "\n"
    "\tporcupine -S pca_s.asc -M 0:3 -T 0.5 model.pdb pca_U.asc >porcupine.pdb\n"
    "Here, a PCA result is assumed, the first 3 modes are used, autoscaling is on,\n"
    "and a 'tip' for the PCA vectors with length 0.5 Angstroms is created.\n"
    "\n"
    "\tporcupine -S pca_s.asc -M 0,3,7 -L 3 model.pdb pca_U.asc >porcupine.pdb\n"
    "A PCA result is assumed, the first, fourth, and eighth mode are used, autoscaling\n"
    "is turned on with a length of 3 Angstroms.  The singular values are also included.\n"
    "\n"
    "\tporcupine --enm -S enm_s.asc -M 0,1 -A 0 --global 50 model.pdb enm_U.asc >porcupine.pdb\n"
    "An ENM result is expected and the first two modes are used.  Autoscaling is disabled.\n"
    "Each mode is scaled by the corresponding eigenvalue (inverted, since this is an ENM).\n"
    "A global scaling of 50 is applied to all modes.\n"
    "\n"
    "SEE ALSO\n"
    "\n"
    "\tPackages/ElasticNetworks/enmovie\n"
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


double averageSubvectorLength(const Matrix& U, const uint col) {
  double avg = 0.0;

  for (uint i=0; i<U.rows(); i += 3) {
    double l = 0.0;
    for (uint j=0; j<3; ++j)
      l += U(i+j, col) * U(i+j, col);
    avg += sqrt(l);
  }

  return(avg / (U.rows()/3));
}


vector<double> determineScaling(const Matrix& U) {
  
  uint n = cols.size();
  vector<double> scaling(n, 1.0);
  vector<double> svals(n, 1.0);
  vector<double> avgs(n, 0.0);

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

  if (autoscale) {
    // Find the largest avg subvector size...
    double maxscale = 0.0;
    for (uint i=0; i<n; ++i) {
      double avg = averageSubvectorLength(U, cols[i]);
      avgs[i] = avg;
      double scale = avg * scaling[i];
      if (scale > maxscale)
        maxscale = scale;
    }

    for (uint i=0; i<n; ++i)
      scaling[i] *= autolength / maxscale;
  }

  // Incorporate additional scaling...
  if (verbosity > 1) {
    cerr << boost::format("%4s %4s %15s %15s %15s\n") % "col" % "mode" % "sval" % "avg" % "scale";
    cerr << boost::format("%4s %4s %15s %15s %15s\n") % "----" % "----" % "---------------" % "---------------" % "---------------";
  }
  for (uint i=0; i<n; ++i) {
    scaling[i] *= scales[i] * global_scale;
    if (verbosity > 1)
      cerr << boost::format("%4d %4d %15.5f %15.5f %15.5f\n") % cols[i] % (cols[i]-offset) % svals[i] % avgs[i] % scaling[i];
    else if (verbosity > 0)
      cerr << boost::format("Scaling column %d by %f\n") % cols[i] % scaling[i];
  }

  return(scaling);
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection* sopts = new opts::BasicSelection("name == 'CA'");
  opts::ModelWithCoords* mopts = new opts::ModelWithCoords;
  ToolOptions* topts = new ToolOptions;
  opts::RequiredArguments *ropts = new opts::RequiredArguments;
  ropts->addArgument("lsv", "left-singular-vector-file");

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mopts).add(topts).add(ropts);
  if (!options.parse(argc, argv))
    exit(-1);

  verbosity = bopts->verbosity;

  // First, read in the LSVs
  Matrix U;
  readAsciiMatrix(ropts->value("lsv"), U);
  uint m = U.rows();

  vector<double> scalings = determineScaling(U);

  // Read in the average structure...
  AtomicGroup avg = mopts->model;

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

  int atomid = 1;
  AtomicGroup spines;

  for (uint j=0; j<cols.size(); ++j) {
    int resid = 1;
    double k = scalings[j];
    uint col = cols[j];

    string segid = generateSegid(col - offset);

    for (uint i=0; i<m; i += 3) {
      GCoord v(U(i,col), U(i+1,col), U(i+2,col));
      if (uniform)
        v /= v.length();

      v *= k;
      pAtom pa = avg.findById(atomids[i/3]);
      GCoord c = pa->coords();

      pAtom atom2;
      if (double_sided)
        atom2 = pAtom(new Atom(atomid++, porcupine_tag, c-v));
      else
        atom2 = pAtom(new Atom(atomid++, porcupine_tag, c));

      atom2->resid(resid);
      atom2->resname(porcupine_tag);
      atom2->segid(segid);


      if (tip_size == 0.0) {
        pAtom atom1(new Atom(atomid++, porcupine_tag, c+v));
        atom1->resid(resid);
        atom1->resname(porcupine_tag);
        atom1->segid(segid);


        atom1->addBond(atom2);
        atom2->addBond(atom1);

        spines.append(atom2);
        spines.append(atom1);

      } else {

        GCoord base = c + v;
        GCoord tip = base + (v / v.length()) * tip_size;

        pAtom atom1(new Atom(atomid++, porcupine_tag, base));
        atom1->resid(resid);
        atom1->resname(porcupine_tag);
        atom1->segid(segid);
        atom1->addBond(atom2);
        atom2->addBond(atom1);

        pAtom atom0(new Atom(atomid++, tip_tag, tip));
        atom0->resid(resid);
        atom0->resname(porcupine_tag);
        atom0->segid(segid);
        atom1->addBond(atom0);
        atom0->addBond(atom1);
        
        spines.append(atom2);
        spines.append(atom1);
        spines.append(atom0);
      }
      ++resid;
    }

  }

  PDB outpdb = PDB::fromAtomicGroup(spines);
  outpdb.remarks().add(hdr);
  cout << outpdb;
}

