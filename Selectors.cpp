/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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

#include <Selectors.hpp>

namespace loos {

  bool CAlphaSelector::operator()(const pAtom& pa) const {
    return(pa->name() == "CA");
  }

  bool BackboneSelector::operator()(const pAtom& pa) const {
    std::string s = pa->name();
    return(s == "C" || s == "CA" || s == "O" || s == "N");
  }

  bool SegidSelector::operator()(const pAtom& pa) const {
    return(pa->segid() == str);
  }

  bool ResidRangeSelector::operator()(const pAtom& pa) const {
    return(pa->resid() >= _low && pa->resid() <= _high);
  }

  
  bool ZSliceSelector::operator()(const pAtom& pa) const {
    greal z = (pa->coords()).z();
    return ( (z>=_min) && (z<_max) );
  }

  bool NotSelector::operator()(const pAtom& pa) const {
    return(!(sel(pa)));
  }

  bool HydrogenSelector::operator()(const pAtom& pa) const {
    bool masscheck = true;

    if (pa->checkProperty(Atom::massbit))
      masscheck = (pa->mass() < 1.1);
    
    std::string n = pa->name();
    return( (n[0] == 'H') && masscheck );
  }

  bool HeavyAtomSelector::operator()(const pAtom& pa) const {
    return (not_heavy(pa));
  }


  bool AndSelector::operator()(const pAtom& pa) const {
    return(lhs(pa) && rhs(pa));
  }

  bool OrSelector::operator()(const pAtom& pa) const {
    return(lhs(pa) || rhs(pa));
  }


  bool SolventSelector::operator()(const pAtom& pa) const {
    return(osel(pa));
  }


  bool HeavySolventSelector::operator()(const pAtom& pa) const {
    return(sel(pa));
  }

  bool KernelSelector::operator()(const pAtom& pa) const {
    krnl.execute(pa);
    if (krnl.stack().size() != 1) {
      throw(std::runtime_error("Execution error - unexpected values on stack"));
    }
    
    internal::Value results = krnl.stack().pop();
    if (results.type != internal::Value::INT)
      throw(std::runtime_error("Execution error - unexpected value on top of stack"));
    
    return(results.itg);
  }

}
