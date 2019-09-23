/*

chain_state_probs.cpp

*/

/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011, Tod D. Romo
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

#include "loos.hpp"

using namespace std;
using namespace loos;

int main(int argc, char *argv[]) {

    string hdr = invocationHeader(argc, argv);
    cout << "# " << hdr << endl;

    char *system_file = argv[1];
    char *traj_file = argv[2];
    char *lipid_selection = argv[3];

    AtomicGroup system = createSystem(system_file);
    pTraj traj = createTrajectory(traj_file, system);

    vector<AtomicGroup> molecules = system.splitByMolecule();
    vector<AtomicGroup> chains;

    for (vector<AtomicGroup>::const_iterator m = molecules.begin();
                                             m != molecules.end();
                                             ++m)
        {
        AtomicGroup tmp = selectAtoms(*m, lipid_selection);
        if (tmp.size() > 0)
            {
            chains.push_back(tmp);
            }
        }

    GCoord normal = GCoord(0.0, 0.0, 1.0);

    uint num_segs = chains[0].size() - 1;

    ChainState states = ChainState(num_segs, 5);

    while (traj->readFrame())
        {
        traj->updateGroupCoords(system);
        for (vector<AtomicGroup>::const_iterator c = chains.begin();
                                                 c != chains.end();
                                                 ++c)
            if (c->centroid().z() > 0.0)
                {
                states.computeChainState(*c, normal);
                }
            else {
                states.computeChainState(*c, -normal);
                }
        }
    std::set<std::pair<StateVector, uint>, Comparator > all_freqs =
            states.getAllProbs();

    cout << "# num_counts = " << states.num_counts() << endl;
    cout << "# Entropy = " << states.entropy() << endl;
    cout << "# Prob\tState" << endl;
    for (std::set<std::pair<StateVector, uint>, Comparator >::iterator
            p = all_freqs.begin();
            p != all_freqs.end();
            ++p)
            {
            double prob = static_cast<double>(p->second) / states.num_counts();
            cout << prob << "\t";
            StateVector s = p->first;
            for (uint i=0; i<s.size(); ++i)
                {
                cout << s[i] << "\t";
                }
            cout << endl;
            }
}
