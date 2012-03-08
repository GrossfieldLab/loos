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
#include <XForm.hpp>

  typedef double                   greal;
  typedef loos::Matrix44<double>   GMatrix;
  typedef loos::Coord<double>      GCoord;

%}


namespace loos {

  typedef Matrix44<greal> GMatrix;

  class XForm {
  public:
    XForm();
    explicit XForm(const GMatrix& m);
    void push(void);
    void pop(void);
    void load(const GMatrix&);
    void concat(const GMatrix&);
    void premult(const GMatrix&);
    void identity(void);
    bool unset(void) const;
    void translate(const greal, const greal, const greal);
    void translate(const GCoord&);
    void scale(const greal, const greal, const greal);
    void scale(const GCoord&);
    void rotate(const GCoord&, const greal);
    void rotate(const char, const greal);
    GCoord transform(const GCoord&);
    GMatrix current(void) const;

  };

  %rename(XFormVector)    std::vector<XForm>;

};
