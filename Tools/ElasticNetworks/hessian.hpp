/*
  hessian

  (c) 2009,2010 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

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


#if !defined(HESSIAN_HPP)
#define HESSIAN_HPP


#include <loos.hpp>

#include "spring_functions.hpp"


//! Base class for defining different methods of adjusting spring constants within the Hessian...

class SuperBlock {
public:
  SuperBlock() : springs(0) { }
  SuperBlock(SpringFunction* func, const loos::AtomicGroup& nodelist) : springs(func), nodes(nodelist) { }
  virtual ~SuperBlock() { }

  uint size() const { return(static_cast<uint>(nodes.size())); }

  virtual loos::DoubleMatrix block(const uint j, const uint i) {
    if (i >= size() || j >= size())
      throw(std::runtime_error("Invalid index in Hessian SuperBlock"));

    if (springs == 0)
      throw(std::runtime_error("No spring function defined for hessian!"));

    loos::GCoord u = nodes[i]->coords();
    loos::GCoord v = nodes[j]->coords();
    loos::GCoord d = v - u;
    
    loos::DoubleMatrix K = springs->constant(u, v, d);
    loos::DoubleMatrix B(3, 3);
    for (uint y=0; y<3; ++y)
      for (uint x=0; x<3; ++x)
        B(x, y) = d[x]*d[y] * K(x,y);

    return(B);
  }


protected:
  SpringFunction* springs;
  loos::AtomicGroup nodes;
};


class SuperBlockDecorator : public SuperBlock {
public:
  SuperBlockDecorator(SuperBlock* b) : decorated(b) { }

protected:
  SuperBlock *decorated;
};



class BoundSuperBlock : public SuperBlockDecorator {
public:
  BoundSuperBlock(SuperBlock* b, SpringFunction* bs, loos::Math::Matrix<int>& cm) :
    SuperBlockDecorator(b),
    bound_spring(bs),
    connectivity(cm) { }

  loos::DoubleMatrix block(const uint j, const uint i) {
    if (j >= connectivity.rows() || i >= connectivity.cols())
      throw(std::runtime_error("Invalid connectivity index"));
    if (j >= size() || i >= size())
      throw(std::runtime_error("Invalid index in Hessian SuperBlock"));
    if (bound_spring == 0)
      throw(std::runtime_error("No spring function defined for bound hessian"));

    loos::DoubleMatrix B(3, 3);
    if (connectivity(j, i)) {
      loos::GCoord u = nodes[i]->coords();
      loos::GCoord v = nodes[j]->coords();
      loos::GCoord d = v - u;

    
      loos::DoubleMatrix K = springs->constant(u, v, d);
      for (uint y=0; y<3; ++y)
        for (uint x=0; x<3; ++x)
          B(x, y) = d[x]*d[y] * K(x,y);
    } else
      B = decorated->block(j, i);

    return(B);
  }


private:
  SpringFunction* bound_spring;
  loos::Math::Matrix<int> connectivity;

};


loos::DoubleMatrix hessian(SuperBlock* block);



#endif
