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


class SuperBlock {
public:
  SuperBlock(const loos::AtomicGroup& nodelist) : nodes(nodelist) { }
  
  loos::DoubleMatrix block(const uint j, const uint i) {
    if (j >= size() || i >= size())
      throw(std::range_error("Invalid indices into SuperBlock"));
    return(blockImpl(j, i));
  }


  uint size() const { return(static_cast<uint>(nodes.size())); }


protected:
  loos::AtomicGroup nodes;

private:
  virtual loos::DoubleMatrix blockImpl(const uint j, const uint i) =0;

};




class DistanceCutoff : public SuperBlock {
public:
  DistanceCutoff(const loos::AtomicGroup& nodelist, const double r) : SuperBlock(nodelist), radius(r*r) { }
private:
  loos::DoubleMatrix blockImpl(const uint j, const uint i);
  double radius;
};


class DistanceWeight : public SuperBlock {
public:
  DistanceWeight(const loos::AtomicGroup& nodelist) : SuperBlock(nodelist), power(-2.0) { }
  DistanceWeight(const loos::AtomicGroup& nodelist, const double power_) : SuperBlock(nodelist), power(power_) { }
private:
  loos::DoubleMatrix blockImpl(const uint j, const uint i);
  double power;
};


class ExponentialDistance : public SuperBlock {
public:
  ExponentialDistance(const loos::AtomicGroup& nodelist) : SuperBlock(nodelist), scale(-1.0) { }
  ExponentialDistance(const loos::AtomicGroup& nodelist, const double scale_) : SuperBlock(nodelist), scale(scale_) { }
private:
  loos::DoubleMatrix blockImpl(const uint j, const uint i);
  double scale;
};


loos::DoubleMatrix hessian(SuperBlock* block);

#endif
