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




%header %{
#include <cryst.hpp>
%}




namespace loos {

  typedef double                greal;

  //! This class encapsulates crystallographic unit cell data
  class UnitCell {
  public:
    UnitCell();
    UnitCell(const GCoord& v);
    greal a(void) const;
    void a(const greal x);

    greal b(void) const;
    void b(const greal x);

    greal c(void) const;
    void c(const greal x);

    greal alpha(void) const;
    void alpha(const greal x);

    greal beta(void) const;
    void beta(const greal x);

    greal gamma(void) const;
    void gamma(const greal x);

    std::string spaceGroup(void) const;
    void spaceGroup(const std::string s);

    int z(void) const;
    void z(const int i);

    %extend {
      char* __str__() {
        std::ostringstream oss;
        oss << *$self;
        size_t n = oss.str().size();
        char* buf = new char[n+1];
        strncpy(buf, oss.str().c_str(), n+1);
        return(buf);
      }

    }

  };


}
