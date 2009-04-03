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





#include <ensembles.hpp>
#include <XForm.hpp>
#include <AtomicGroup.hpp>
#include <Trajectory.hpp>

namespace loos {

  // Assume all groups are already sorted or matched...

  AtomicGroup averageStructure(const std::vector<AtomicGroup>& ensemble) {
    AtomicGroup avg = ensemble[0].copy();

    // First, zap our coords...
    int n = avg.size();
    int i;
    for (i=0; i<n; i++)
      avg[i]->coords() = GCoord(0.0, 0.0, 0.0);

    // Now, accumulate...
    std::vector<AtomicGroup>::const_iterator j;
    for (j = ensemble.begin(); j != ensemble.end(); j++) {
      for (i = 0; i<n; i++)
        avg[i]->coords() += (*j)[i]->coords();
    }

    for (i=0; i<n; i++)
      avg[i]->coords() /= ensemble.size();

    return(avg);
  }



  AtomicGroup averageStructure(const AtomicGroup& g, const std::vector<XForm>& xforms, pTraj traj) {
    AtomicGroup avg = g.copy();
    AtomicGroup frame = g.copy();
    int n = avg.size();
    for (int i=0; i<n; i++)
      avg[i]->coords() = GCoord(0.0, 0.0, 0.0);

    traj->rewind();
    uint tn = traj->nframes();
    
    // Safety check...
    if (tn != xforms.size())
      throw(std::runtime_error("Mismatch in number of frames in the trajectory and passed transforms for loos::averageStructure()"));

    for (uint j=0; j<tn; j++) {
      traj->readFrame(j);
      traj->updateGroupCoords(frame);
      frame.applyTransform(xforms[j]);
      for (int i=0; i<n; i++)
        avg[i]->coords() += frame[i]->coords();
    }

    for (int i=0; i<n; i++)
      avg[i]->coords() /= tn;

    return(avg);
  }




  boost::tuple<std::vector<XForm>,greal,int> iterativeAlignment(std::vector<AtomicGroup>& ensemble,
                                                                      greal threshold, int maxiter) {
    int iter = 0;
    int n = ensemble.size();
    greal rms;
    std::vector<XForm> xforms(n);
    AtomicGroup avg;

    // Start by aligning against the first structure in the ensemble
    AtomicGroup target = ensemble[0].copy();

    target.centerAtOrigin();

    do {
      for (int i = 0; i<n; i++) {
        GMatrix M = ensemble[i].alignOnto(target);
        xforms[i].premult(M);
      }

      avg = averageStructure(ensemble);
      rms = avg.rmsd(target);
      target = avg;
      iter++;
    } while (rms > threshold && iter <= maxiter );
  
    boost::tuple<std::vector<XForm>, greal, int> res(xforms, rms, iter);
    return(res);
  }




  boost::tuple<std::vector<XForm>, greal, int> iterativeAlignment(const AtomicGroup& g, pTraj traj,
                                                                        greal threshold, int maxiter) {
    std::vector<AtomicGroup> frames;

    traj->rewind();
    while (traj->readFrame()) {
      AtomicGroup frame = g.copy();
      traj->updateGroupCoords(frame);
      frames.push_back(frame);
    }

    boost::tuple<std::vector<XForm>, greal, int> result = iterativeAlignment(frames, threshold, maxiter);
    return(result);
  }


  void applyTransforms(std::vector<AtomicGroup>& ensemble, std::vector<XForm>& xforms) {
    uint n = ensemble.size();
    if (n != xforms.size())
      throw(std::runtime_error("Mismatch in the size of the ensemble and the transformations"));

    for (uint i=0; i<n; ++i)
      ensemble[i].applyTransform(xforms[i]);
  }



  void readTrajectory(std::vector<AtomicGroup>& ensemble, const AtomicGroup& model, pTraj trajectory) {
    AtomicGroup clone = model.copy();
    
    while (trajectory->readFrame()) {
      trajectory->updateGroupCoords(clone);
      AtomicGroup frame = clone.copy();
      ensemble.push_back(frame);
    }
  }


  void readTrajectory(std::vector<AtomicGroup>& ensemble, const AtomicGroup& model, pTraj trajectory, std::vector<uint>& frames) {
    AtomicGroup clone = model.copy();
    
    std::vector<uint>::iterator i;
    for (i = frames.begin(); i != frames.end(); ++i) {
      if (*i >= trajectory->nframes())
        throw(std::runtime_error("Frame index exceeds trajectory size in readTrajectory()"));
      trajectory->readFrame(*i);
      trajectory->updateGroupCoords(clone);
      AtomicGroup frame = clone.copy();
      ensemble.push_back(frame);
    }
  }



  RealMatrix extractCoords(std::vector<AtomicGroup>& ensemble) {
    uint n = ensemble.size();
    uint m = ensemble[0].size() * 3;
    RealMatrix M(m, n);

    for (uint i=0; i<n; ++i)
      for (uint j=0; j<m; ++j) {
        M(3*j, i) = ensemble[i][j]->coords().x();
        M(3*j+1, i) = ensemble[i][j]->coords().y();
        M(3*j+2, i) = ensemble[i][j]->coords().z();
      }

    return(M);
  }


  RealMatrix extractCoords(std::vector<AtomicGroup>& ensemble, std::vector<XForm>& xforms) {
    uint n = ensemble.size();

    if (n != xforms.size())
      throw(std::runtime_error("Mismatch between the size of the ensemble and the transformations"));

    uint m = ensemble[0].size() * 3;
    RealMatrix M(m, n);

    for (uint i=0; i<n; ++i) {
      GMatrix W = xforms[i].current();

      for (uint j=0; j<m; ++j) {
        GCoord c = W * ensemble[i][j]->coords();
        M(3*j, i) = c.x();
        M(3*j+1, i) = c.y();
        M(3*j+2, i) = c.z();
      }
    }

    return(M);
  }


  void subtractAverage(RealMatrix& M) {
    uint m = M.rows();
    uint n = M.cols();
    double *avg = new double[m];
    
    for (uint j=0; j<m; ++j)
      avg[j] = 0.0;

    for (uint i=0; i<n; ++i)
      for (uint j=0; j<m; ++j)
        avg[j] += M(j, i);

    for (uint j=0; j<m; ++j)
      avg[j] /= n;

    for (uint i=0; i<n; ++i)
      for (uint j=0; j<m; ++j)
        M(j,i) -= avg[j];
  }

  boost::tuple<RealMatrix, RealMatrix, RealMatrix> svd(RealMatrix& M) {
    f77int m = M.rows();
    f77int n = M.cols();
    f77int sn = m<n ? m : n;
    long estimate = (m*m + n*n + m*n + sn) * sizeof(float);
    // This is somewhat vestigial...  We continue to estimate the
    // total storage required for the SVD in case we want to do
    // something with it in the future...



    char jobu = 'A', jobvt = 'A';
    f77int lda = m, ldu = m, ldvt = n, lwork= -1, info;
    float prework[10], *work;
    
    RealMatrix U(m, m);
    RealMatrix S(sn, 1);
    RealMatrix Vt(n, n);

    sgesvd_(&jobu, &jobvt, &m, &n, M.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, prework, &lwork, &info);
    if (info != 0) {
      // throw here...
      exit(-2);
    }

    lwork = (f77int)prework[0];
    estimate += lwork * sizeof(double);
    work = new float[lwork];

    sgesvd_(&jobu, &jobvt, &m, &n, M.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, work, &lwork, &info);

    boost::tuple<RealMatrix, RealMatrix, RealMatrix> result(U, S, Vt);
    return(result);
  }



  boost::tuple<DoubleMatrix, DoubleMatrix, DoubleMatrix> svd(DoubleMatrix& M) {
    f77int m = M.rows();
    f77int n = M.cols();
    f77int sn = m<n ? m : n;
    long estimate = (m*m + n*n + m*n + sn) * sizeof(double);


    char jobu = 'A', jobvt = 'A';
    f77int lda = m, ldu = m, ldvt = n, lwork= -1, info;
    double prework[10], *work;
    
    DoubleMatrix U(m, m);
    DoubleMatrix S(sn, 1);
    DoubleMatrix Vt(n, n);

    dgesvd_(&jobu, &jobvt, &m, &n, M.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, prework, &lwork, &info);
    if (info != 0) {
      // throw here...
      exit(-2);
    }

    lwork = (f77int)prework[0];
    estimate += lwork * sizeof(double);
    work = new double[lwork];

    dgesvd_(&jobu, &jobvt, &m, &n, M.get(), &lda, S.get(), U.get(), &ldu, Vt.get(), &ldvt, work, &lwork, &info);

    boost::tuple<DoubleMatrix, DoubleMatrix, DoubleMatrix> result(U, S, Vt);
    return(result);
  }


  boost::tuple<RealMatrix, RealMatrix, RealMatrix> svd(std::vector<AtomicGroup>& ensemble, bool align = true) {

    RealMatrix M;
    if (align) {
      boost::tuple<std::vector<XForm>, greal, int> res = iterativeAlignment(ensemble);
      M = extractCoords(ensemble, boost::get<0>(res));
    } else
      M = extractCoords(ensemble);

    subtractAverage(M);
    boost::tuple<RealMatrix, RealMatrix, RealMatrix> res = svd(M);
    return(res);
  }

}
