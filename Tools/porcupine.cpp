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
vector<uint> cols;
vector<double> scales;
double global_scale;
bool uniform;
bool double_sided;
string model_name;
string map_name;
double tip_size;
string selection;

const string porcupine_tag("POR");
const string tip_tag("POT");

// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("columns,C", po::value< vector<string> >(&strings), "Columns to use")
      ("scale,S", po::value< vector<double> >(&scales), "Scale the requested columns")
      ("global", po::value<double>(&global_scale)->default_value(1.0), "Global scaling")
      ("uniform", po::value<bool>(&uniform)->default_value(false), "Scale all elements uniformly")
      ("map", po::value<string>(&map_name), "Use a map file to map LSV/eigenvectors to atomids")
      ("tips,T", po::value<double>(&tip_size)->default_value(0.0), "Length (in Angstroms) to make the tip (for single-sided only)")
      ("double_sided", po::value<bool>(&double_sided)->default_value(false), "Use double-sided vectors");
  }

  bool postConditions(po::variables_map& vm) {
    if (strings.empty())
      cols.push_back(0);
    else
      cols = parseRangeList<uint>(strings);

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
    oss << boost::format("columns='%s', global=%f, uniform=%d, map='%s', tips=%f, double_sided=%d")
      % vectorAsStringWithCommas<string>(strings)
      % global_scale
      % uniform
      % map_name
      % tip_size
      % double_sided;
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
    "EXAMPLES\n"
    "\n"
    "\tporcupine -s 'name==\"CA\"' -C1 -S2 --double_sided 1 model.pdb svd_U.asc > porky.pdb\n"
    "Here the CA's from model.pdb are used as the origin of the vector drawing.\n"
    "svd_U.asc is an LSV file made in loos (see example X in svd --fullhelp)\n"
    "These vectors are mapped to the CA's in model.pdb.  The C1 option signifies\n"
    "use of the 1st column...in this case the most collective motion from the PCA.\n"
    "The option -S2 implies an arbitrary scaling factor of 2 is applied to the LSVs\n"
    "Double-sided vectors are turned on - this means the vectors will be drawn in\n"
    "both directions from the CA origin.  Finally, a new pdb, porky.pdb is created.\n"
    "This pdb contains only the vectors drawn by the porcupine program.\n"
    "\t---You may wish to visualize the output of this tool with pymol.\n"
    "\t   Try using the following options (within pymol):\n"
    "\t\t load model.pdb; load porky.pdb\n"
    "\t\t hide\n"
    "\t\t show cartoon, model\n"
    "\t\t show sticks, porky\n"
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



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* basopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::ModelWithCoords* mopts = new opts::ModelWithCoords;
  ToolOptions* topts = new ToolOptions;
  opts::RequiredArguments *ropts = new opts::RequiredArguments;
  ropts->addArgument("lsv", "left-singular-vector-file");

  opts::AggregateOptions options;
  options.add(basopts).add(bopts).add(sopts).add(mopts).add(topts).add(ropts);
  if (!options.parse(argc, argv))
    exit(-1);

  // First, read in the LSVs
  Matrix U;
  readAsciiMatrix(ropts->value("lsv"), U);
  uint m = U.rows();

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
  int resid = 1;
  AtomicGroup spines;

  for (uint j=0; j<cols.size(); ++j) {
    double k = global_scale * scales[j];
    uint col = cols[j];

    string segid = generateSegid(col);

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

      atom2->resid(resid++);
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
    }

  }

  PDB outpdb = PDB::fromAtomicGroup(spines);
  outpdb.remarks().add(hdr);
  cout << outpdb;
}

