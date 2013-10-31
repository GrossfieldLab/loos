/*
  phase_pdb

  Takes 3 colums from a matrix and sticks them into the coordinates for a fake PDB.
  This is an easy way to visualize PCA results in 3D using your
  favorite molecular visualization software...
*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010, Tod D. Romo
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


vector<uint> columns;
vector<uint> rows;
vector<double> scales;
uint chunksize;
string segid_fmt;
string residue_name;
string atom_name;

bool bonds;

string matrix_name;


// @cond TOOLS_INTERNAL


string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tCreates pseudoatoms representing the phase-space of a trajectory\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tphase-pdb takes right singular vectors (RSV) from an SVD of a trajectory and creates\n"
    "pseudoatoms at each point in phase-space.  This PDB can be loaded into pymol or vmd\n"
    "and used to visualize the phase-space projection of the trajectory.  In order to\n"
    "show the time-evolution of the system, the pseudoatoms can be connected by bonds using\n"
    "the --bonds=1 option.\n"
    "\tThe elements of the RSV will need to be scaled up in order to be visualized.  This is\n"
    "done with the --scales option.  Additionally, to visualize the true shape of the phase-\n"
    "space, the RSV columns should be scaled by the corresponding singular values.  This must\n"
    "be done manually, i.e. find the appropriate singular values, scale them by a constant,\n"
    "and use this with the --scales option.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tphase-pdb b2ar_V.asc >b2ar_V.pdb\n"
    "This uses the first 3 RSVs, scaling each with a default of 100.\n"
    "\n"
    "\tphase-pdb --scales=100 --scales=50 --scales=25 b2ar_V.asc >b2ar_V.pdb\n"
    "This uses the first 3 RSVs, scaling the first by 100, the second by 50, and the third\n"
    "by 25.\n"
    "\n"
    "\tphase-pdb --bonds=1 --scales=100 --scales=50 --scales=25 b2ar_V.asc >b2ar_V.pdb\n"
    "This uses the first 3 RSVs, scaling them as above, but adds CONECT records at the\n"
    "end of the PDB that connects sequential pseudoatoms.\n"
    "\n"
    "SEE ALSO\n"
    "\tsvd, big-svd\n";

  return(msg);
}


class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("segid", po::value<string>(&segid_fmt)->default_value("P%03d"), "Segid printf format")
      ("atom", po::value<string>(&atom_name)->default_value("CA"), "Atom name to use")
      ("residue", po::value<string>(&residue_name)->default_value("SVD"), "Residue name to use")
      ("rows", po::value<string>(&rowdesc)->default_value("all"), "Rows to use")
      ("column,C", po::value< vector<uint> >(&columns), "Columns to use (default are first 3)")
      ("scales,S", po::value< vector<double> >(&scales), "Scale columns (default is 100,100,100)")
      ("chunk", po::value<uint>(&chunksize)->default_value(0), "Divide vector into chunks by these number of rows")
      ("bonds", po::value<bool>(&bonds)->default_value(false), "Connect sequential atoms by bonds");
    
  }

  bool postConditions(po::variables_map& vm) {
    if (columns.empty())
      for (uint i=0; i<3; ++i)
        columns.push_back(i);
    
    if (scales.empty())
      for (uint i=0; i<3; ++i)
        scales.push_back(100.0);
    else if (scales.size() == 1)
      for (uint i=0; i<2; ++i)
        scales.push_back(scales[0]);

    if (columns.size() != scales.size()) {
      cerr << "Error- number of columns selected does not equal number of scales\n";
      return(false);
    }

    if (columns.size() != 3) {
      cerr << "Error- must select 3 columns\n";
      return(false);
    }

    if (rowdesc != "all")
      rows = parseRangeList<uint>(rowdesc);
    
    return(true);
  }


  string print() const {
    ostringstream oss;
    oss << boost::format("segid='%s', atom='%s', residue='%s', rows='%s', chunk=%d, bonds=%d, columns=(%s), scales=(%f)")
      % segid_fmt
      % atom_name
      % residue_name
      % rowdesc
      % rowdesc
      % chunksize
      % bonds
      % vectorAsStringWithCommas<uint>(columns)
      % vectorAsStringWithCommas<double>(scales);

    return(oss.str());
  }

  string rowdesc;

};

// @endcond




int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  ToolOptions* topts = new ToolOptions;
  opts::RequiredArguments* ropts = new opts::RequiredArguments;
  ropts->addArgument("matrix", "matrix-file");

  opts::AggregateOptions options;
  options.add(bopts).add(topts).add(ropts);
  if (!options.parse(argc, argv))
    exit(-1);

  string matrix_name = ropts->value("matrix");

  RealMatrix A;
  readAsciiMatrix(matrix_name, A);
  uint m = A.rows();

  if (rows.empty())
    for (uint i=0; i<m; ++i)
      rows.push_back(i);

  uint resid = 1;
  uint total_chunks = rows.size() / chunksize;
  uint chunk = 1;
  AtomicGroup model;

  for (uint atomid = 0; atomid < rows.size(); ++atomid, ++resid) {
    if (chunksize && resid > chunksize) {
      if (bonds) {
        for (uint i=atomid - chunksize; i<atomid-1; ++i)
          model[i]->addBond(model[i+1]);
      }

      resid = 1;
      ++chunk;
    }

    GCoord c(scales[0] * A(rows[atomid], columns[0]),
             scales[1] * A(rows[atomid], columns[1]),
             scales[2] * A(rows[atomid], columns[2]));
    pAtom pa(new Atom(atomid+1, atom_name, c));
    pa->resid(resid);
    pa->resname(residue_name);

    ostringstream segstream;
    segstream << boost::format(segid_fmt) % chunk;
    pa->segid(segstream.str());

    double b;
    double q = 0.0;
    
    if (chunksize) {
      b = (100.0 * (resid-1)) / chunksize;
      q = static_cast<double>(chunk-1) / total_chunks;
    } else
      b = 100.0 * atomid / rows.size();

    pa->bfactor(b);
    pa->occupancy(q);

    model.append(pa);
  }
  
  if (bonds) {
      if (chunksize && resid > 1)
	  for (uint i=rows.size() - chunksize; i<rows.size()-1; ++i)
	      model[i]->addBond(model[i+1]);
      else if (!chunksize)
	  for (uint i=0; i<rows.size()-1; ++i)
	      model[i]->addBond(model[i+1]);
  }
  
  PDB pdb = PDB::fromAtomicGroup(model);
  pdb.remarks().add(hdr);
  cout << pdb;
}
