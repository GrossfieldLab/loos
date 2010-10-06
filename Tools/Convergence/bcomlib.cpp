/*
  bcomlib.cpp

  (c) 2010 Tod D. Romo, Grossfield Lab, URMC
*/



#include "bcomlib.hpp"



void Convergence::subtractStructure(loos::RealMatrix& M, const loos::AtomicGroup& model) {
  std::vector<float> avg(model.size() * 3);
  int k = 0;
  for (int i=0; i<model.size(); ++i) {
    loos::GCoord c = model[i]->coords();
    avg[k++] = c.x();
    avg[k++] = c.y();
    avg[k++] = c.z();
  }

  for (uint i=0; i<M.cols(); ++i)
    for (uint j=0; j<M.rows(); ++j)
      M(j, i) -= avg[j];
}
