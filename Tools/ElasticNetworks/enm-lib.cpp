/*
  enm-lib

  (c) 2009 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry


  Common code for the ENM suite

*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009 Tod D. Romo
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


#include <enm-lib.hpp>


using namespace std;
using namespace loos;


AtomicGroup sideChainCentroids(const AtomicGroup& grp, int maxid, int maxresid, const string& name, const string& resname, const string& segid) {
  
  if (maxid == 0)
    maxid = grp.maxId() + 1;
  if (maxresid == 0)
    maxresid = grp.maxResid() + 1;

  // Build the non-backbone selector
  AtomNameSelector csel("C");
  AtomNameSelector casel("CA");
  OrSelector or1(csel,casel);

  AtomNameSelector nsel("N");
  OrSelector or2(or1, nsel);

  AtomNameSelector osel("O");
  OrSelector or3(or2, osel);

  HydrogenSelector hsel;
  OrSelector or4(or3, hsel);

  NotSelector notbackbone(or4);

  AtomicGroup result;
  vector<AtomicGroup> residues = grp.splitByResidue();
  for (AtomicGroup::iterator i = residues.begin(); i != residues.end(); ++i) {
    AtomicGroupy sidechain = (*i)->select(notbackbone);
    if (sidechain.empty())
      continue;

    GCoord center = sidechain.centroid();
    greal mass = sidechain.totalMass();

    pAtom pa(new Atom(maxid++, name, center));
    pa->resid(maxresid++);
    pa->resname(resname);
    pa->segid(segid);
    pa->mass(mass);
    result.append(pa);
  }

  return(result);
}
