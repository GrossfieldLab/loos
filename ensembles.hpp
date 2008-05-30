/*
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Ensemble calculations...
*/

#if !defined(ENSEMBLES_HPP)
#define ENSEMBLES_HPP

namespace loos {
  AtomicGroup averageStructure(const vector<AtomicGroup>&);
  greal iterativeAlignment(vector<AtomicGroup>& ensmble, greal threshold, int maxiter=1000);
};



#endif
