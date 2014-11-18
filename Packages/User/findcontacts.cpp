// (c) 2014 Tod D. Romo, Grossfield Lab, URMC

#include "loos.hpp"

using namespace std;
using namespace loos;

typedef vector<AtomicGroup>    vGroup;

const double prunefactor = 2.0;



vector<uint> findContacts(const AtomicGroup& source, const GCoord& source_center, const AtomicGroup& target, const double cutoff2, const double prunefactor2) 
{
  vector<uint> contacts(source.size(), 0);

  for (AtomicGroup::const_iterator cj = target.begin(); cj != target.end(); ++cj) {
    GCoord x = (*cj)->coords();
    if (x.distance2(source_center) >= prunefactor2)
      continue;
    for (uint i=0; i<source.size(); ++i)
      if (!contacts[i] && x.distance2(source[i]->coords()) <= cutoff2)
	contacts[i] = 1;
  }

  return(contacts);
}



int main(int argc, char *argv[]) 
{
  if (argc < 8) {
    cerr << "Usage- " << argv[0] << " cutoff model retinal protein water salt traj [traj ...]\n";
    cerr << "Example:\n";
    cerr << "  findcontacts 2.5 npgt_start.pdb 'segid == \"RTNE\" && (hydrogen || name == \"NZ\")' 'segid == \"RHOD\"' 'segid == \"BULK\"' 'segid == \"CHLO\" || segid == \"SODI\"' sim2_1ns.dcd >foo.asc\n";
    exit(-1);
  }
  
  
  
  string hdr = invocationHeader(argc, argv);

  int k = 1;
  double cutoff = strtod(argv[k++], 0);
  AtomicGroup model = createSystem(argv[k++]);
  AtomicGroup retinal = selectAtoms(model, argv[k++]);
  AtomicGroup prot = selectAtoms(model, argv[k++]);
  AtomicGroup water = selectAtoms(model, argv[k++]);
  AtomicGroup salt = selectAtoms(model, argv[k++]);
  
  vGroup residues = prot.splitByResidue();
  vector<string> coltags;
  vector<string> rowtags;

  for (vGroup::iterator i = residues.begin(); i != residues.end(); ++i) {
    ostringstream oss;
    oss << (*i)[0]->segid() << ':' << (*i)[0]->resid();
    coltags.push_back(oss.str());
  }
  
  residues.push_back(water);
  coltags.push_back("water");
  
  residues.push_back(salt);
  coltags.push_back("salt");
  
  for (AtomicGroup::iterator i = retinal.begin(); i != retinal.end(); ++i)
    rowtags.push_back((*i)->name());


  uint n = residues.size();
  uint m = retinal.size();

  vector<uint> colmask(n, 0);

  vector< vector<uint> > C(m, vector<uint>(n, 0));

  cutoff *= cutoff;

  uint t = 0;
  cerr << "Working- ";

  while (k < argc) {
    pTraj traj = createTrajectory(argv[k++], model);
    
    while (traj->readFrame()) {
      if (t++ % 500 == 0)
	cerr << '.';

      traj->updateGroupCoords(model);
      GCoord c = retinal.centroid();
      double r = retinal.radius();
      double prune = r * prunefactor;
      prune *= prune;
    
      for (uint i=0; i<n; ++i) {
	vector<uint> contacts = findContacts(retinal, c, residues[i], cutoff, prune);
	for (uint j=0; j<m; ++j)
	  if (contacts[j]) {
	    colmask[i] = 1;
	    C[j][i] += 1;
	  }
      }
    }
  }
  
  cout << "# " << hdr << endl;
  for (uint j=0; j<m; ++j) {
    uint col = 0;
    
    for (uint i=0; i<n; ++i) {
      if (colmask[i]) {
	cout << j << '\t' << col++ << '\t' << rowtags[j] << '\t' << coltags[i] << '\t' << 
	  (static_cast<double>(C[j][i])/t) << endl;
//	cout << j << '\t' << col++ << '\t' << rowtags[j] << '\t' << coltags[i] << '\t' << 
//	  C[j][i] << endl;

      }
    }
    cout << endl;
  }
  
  cerr << " done" << endl;
  
}
  
