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
#include <Matrix44.hpp>
%}





namespace loos {

  // Forward declaration for matrix-vector multiply
  template<class T> Coord<T> operator*(const Matrix44<T>&, const Coord<T>&);


  //! Specialized 4x4 Matrix class for handling coordinate transforms.
  template<class T>
  class Matrix44 {
  public:
    Matrix44();
    explicit Matrix44(const T v);
    void zero(void);
    void identity(void);
    T& operator()(const int j, const int i);
    T* data(void);
    Matrix44<T>& operator+=(const Matrix44<T>& rhs);
    Matrix44<T> operator+(const Matrix44<T>& rhs);
    Matrix44<T>& operator-=(const Matrix44<T>& rhs);
    Matrix44<T> operator-(const Matrix44<T>& rhs);
    Matrix44<T>& operator*=(const Matrix44<T>& rhs);
    Matrix44<T> operator*(const Matrix44<T>& rhs) const;
    Matrix44<T>& operator*=(const T x);
    Matrix44<T> operator*(const T x);



  };


  %extend Matrix44<double> {
    double __getitem__(const int i) {
      return((*$self)[i]);
    }

    void __setitem__(const int i, const double d) {
      (*$self)[i] = d;
    }

  };

};




%template(GMatrix)   loos::Matrix44<double>;

