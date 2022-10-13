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





#if !defined(LOOS_CRYST_HPP)
#define LOOS_CRYST_HPP

#include <loos_defs.hpp>
#include <Coord.hpp>

namespace loos {

  //! This class encapsulates crystallographic unit cell data
  class UnitCell {
  public:
    UnitCell() : _a(1.0), _b(1.0), _c(1.0),
                 _alpha(90.0), _beta(90.0), _gamma(90.0),
                 sgroup("P1"), zval(1) { }

    UnitCell(const GCoord& v) : _a(v.x()), _b(v.y()), _c(v.z()),
                                _alpha(90.0), _beta(90.0), _gamma(90.0),
                                sgroup("P1"), zval(1) { }

    greal a(void) const { return(_a); }
    void a(const greal x) { _a = x; }

    greal b(void) const { return(_b); }
    void b(const greal x) { _b = x; }

    greal c(void) const { return(_c); }
    void c(const greal x) { _c = x; }

    greal alpha(void) const { return(_alpha); }
    void alpha(const greal x) { _alpha = x; }

    greal beta(void) const { return(_beta); }
    void beta(const greal x) { _beta = x; }

    greal gamma(void) const { return(_gamma); }
    void gamma(const greal x) { _gamma = x; }

    std::string spaceGroup(void) const { return(sgroup); }
    void spaceGroup(const std::string s) { sgroup = s; }

    int z(void) const { return(zval); }
    void z(const int i) { zval = i; }

    friend std::ostream& operator<<(std::ostream& os, const UnitCell& u) {
      os << "<UNITCELL A='" << u._a << "' B='" << u._b << "' C='" << u._c << "' ALPHA='";
      os << u._alpha << "' BETA='" << u._beta << "' GAMMA='" << u._gamma << "' SPACEGROUP='";
      os << u.sgroup << "' Z='" << u.zval << "'/>";
      return(os);
    }


  private:
    greal _a, _b, _c, _alpha, _beta, _gamma;
    std::string sgroup;
    int zval;
  };


}
#endif
