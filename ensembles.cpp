/*
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Ensemble calculations...
*/



#include <loos.hpp>
#include <AtomicGroup.hpp>
#include <ensembles.hpp>


// Assume all groups are already sorted or matched...

AtomicGroup loos::averageStructure(const vector<AtomicGroup>& ensemble) {
  AtomicGroup avg = *(ensemble[0].clone());

  // First, zap our coords...
  int n = avg.size();
  int i;
  for (i=0; i<n; i++)
    avg[i]->coords() = GCoord(0.0, 0.0, 0.0);

  // Now, accumulate...
  vector<AtomicGroup>::const_iterator j;
  for (j = ensemble.begin(); j != ensemble.end(); j++) {
    for (i = 0; i<n; i++)
      avg[i]->coords() += j->getAtomsTransformedCoord(i);
  }

  for (i=0; i<n; i++)
    avg[i]->coords() /= ensemble.size();

  return(avg);
}



greal loos::iterativeAlignment(vector<AtomicGroup>& ensemble, greal threshold, int maxiter) {
  int iter = 0;
  int n = ensemble.size();
  greal rms;
  AtomicGroup avg;
  AtomicGroup target = averageStructure(ensemble);

  do {
    for (int i = 0; i<n; i++)
      ensemble[i].alignOnto(target);
    avg = averageStructure(ensemble);
    rms = avg.rmsd(target);
    target = avg;

#if defined(DEBUG)
    cerr << "loos::iterativeAlignment - iter = " << iter << ", rms = " << rms << endl;
#endif

  } while (rms > threshold && ++iter <= maxiter);
  
  return(rms);
}
