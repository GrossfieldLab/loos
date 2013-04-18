// (c) 2013 Tod D. Romo, Grossfield Lab, URMC


#include <loos.hpp>

using namespace std;
using namespace loos;


int main(int argc, char *argv[]) 
{
    if (argc != 5) {
        cerr << "Usage- serialize-selection output-prefix model trajectory selection\n";
        exit(-1);
    }

    string hdr = invocationHeader(argc, argv);
    
    int k = 1;
    string prefix(argv[k++]);
    AtomicGroup model = createSystem(argv[k++]);
    pTraj traj = createTrajectory(argv[k++], model);
    AtomicGroup subset = selectAtoms(model, argv[k++]);

    vector<AtomicGroup> molecules;
    if (model.hasBonds())
        molecules = subset.splitByMolecule();
    else
        molecules = subset.splitByUniqueSegid();

    bool structure_written = false;
    AtomicGroup outgroup = molecules[0].copy();
    DCDWriter dcd(prefix + ".dcd");

    while (traj->readFrame()) {
        traj->updateGroupCoords(model);

        for (uint j=0; j<molecules.size(); ++j) {
            for (uint i=0; i<outgroup.size(); ++i)
                outgroup[i]->coords(molecules[j][i]->coords());
            if (!structure_written) {
                PDB pdb = PDB::fromAtomicGroup(outgroup);
                pdb.remarks().add(hdr);

                string pdbname = prefix + ".pdb";
                ofstream ofs(pdbname.c_str());
                ofs << pdb;

                structure_written = true;
            }

            dcd.writeFrame(outgroup);
        }
    }
}
