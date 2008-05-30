#include <loos.hpp>
#include <pdb.hpp>
#include <ensembles.hpp>
#include <Selectors.hpp>

const int maxf = 4;


int main() {
  vector<PDB> pdbs;
  int i;

  for (i=0; i<maxf; i++) {
    char buf[256];
    sprintf(buf, "frame_%02d.pdb", i);
    PDB pdb(buf);
    pdbs.push_back(pdb);
  }

  CAlphaSelector casel;
  vector<AtomicGroup> backbones;
  for (i=0; i<maxf; i++) {
    AtomicGroup grp = pdbs[i].select(casel);
    backbones.push_back(grp);
  }

  greal rms = loos::iterativeAlignment(backbones, 0.2);
  cout << "rms = " << rms << endl;

  for (i=0; i<maxf; i++) {
    char buf[256];
    sprintf(buf, "A-frame_%02d.pdb", i);
    pdbs[i].xform().load(backbones[i].xform().current());
    pdbs[i].applyTransform();
    
    ofstream output(buf);
    output << pdbs[i];

  }

}
