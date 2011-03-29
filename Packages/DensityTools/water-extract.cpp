/*
  water-extract.cpp

  (c) 2008 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry


   usage:
     water-extract [options] pdb trajectory matrix >output.pdb

   Desc:
   Given a water matrix and a trajectory/model, will extract the
   internal water atoms, concatenating them all into one big PDB
   file.  This can be used as a quick and dirty way to visualize the
   accumulated water distribution/densities...

*/

#include <iterator>
#include <loos.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>



using namespace std;
using namespace loos;
namespace po = boost::program_options;

typedef Math::Matrix<int, Math::RowMajor>   Matrix;


string water_string, pdb_name, trajectory_name, matrix_name;
vector<uint> cols;


void parseArgs(int argc, char *argv[]) {
  string range;

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("water,w", po::value<string>(&water_string)->default_value("name == 'OH2'"), "Water selection")
      ("range,r", po::value<string>(&range), "Range of columns to extract");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("pdb", po::value<string>(&pdb_name), "Model filename")
      ("traj", po::value<string>(&trajectory_name), "Trajectory filename")
      ("mat", po::value<string>(&matrix_name), "Matrix filename");
    

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("pdb", 1);
    p.add("traj", 1);
    p.add("mat", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("pdb") && vm.count("traj") && vm.count("mat"))) {
      cerr << "Usage- water-extract [options] pdb trajectory matrix >output\n";
      cerr << generic;
      exit(-1);
    }

    if (vm.count("range"))
      cols = parseRange<uint>(range);
        
  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}




int main(int argc, char *argv[])
{
  string hdr = invocationHeader(argc, argv);
  parseArgs(argc, argv);

  AtomicGroup model = createSystem(pdb_name);
  pTraj traj = createTrajectory(trajectory_name, model);

  AtomicGroup water = selectAtoms(model, water_string);

  Matrix M;
  readAsciiMatrix(matrix_name, M);
  uint m = M.rows();
  uint n = M.cols();
  cerr << boost::format("Read in %u x %u matrix from %s.\n") % m % n % matrix_name;

  if (n != traj->nframes()) {
    cerr << boost::format("ERROR - trajectory has %d frames, but expected %d.\n") % traj->nframes() % n;
    exit(-1);
  }
  if (m != water.size()) {
    cerr << boost::format("ERROR - selection results in %d waters, but expected %d.\n") % water.size() % m;
    exit(-1);
  }

  if (cols.empty())
    for (uint i=0; i<n; ++i)
      cols.push_back(i);

  AtomicGroup liquid;
  uint currid = 1;

  for (uint i=0; i<cols.size(); ++i) {
    traj->readFrame(cols[i]);
    traj->updateGroupCoords(water);

    for (uint j=0; j<m; ++j) {
      if (M(j, i)) {
        pAtom atom(new Atom(*(water[j])));
        atom->id(currid);
        atom->resid(currid);
        ++currid;
        liquid.append(atom);
      }
    }

  }
  cerr << boost::format("Concatenated a total of %d waters.\n") % (currid-1);
  PDB pdb = PDB::fromAtomicGroup(liquid);
  pdb.remarks().add(hdr);
  cout << pdb;

}
