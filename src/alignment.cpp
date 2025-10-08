#include <loos_defs.hpp>

#include <AtomicGroup.hpp>
#include <Trajectory.hpp>

#include <ensembles.hpp>
#include <alignment.hpp>

#include <cmath>


namespace loos {



  namespace alignment {


    // Core aligmnent routine.  Assumes input coord vectors are already centered.
    // Returns the SVD results as a tuple
    SVDTupleVec kabschCore(const vecDouble& u, const vecDouble& v) {
      int n = u.size() / 3;

      // Compute correlation matrix...
      vecDouble R(9);

#if defined(__linux__) || defined(__CYGWIN__) || defined(__FreeBSD__)
      char ta = 'N';
      char tb = 'T';
      f77int three = 3;
      double one = 1.0;
      double zero = 0.0;

      dgemm_(&ta, &tb, &three, &three, &n, &one, u.data(), &three, v.data(), &three, &zero, R.data(), &three);

#else

      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 3, 3, n, 1.0, u.data(), 3, v.data(), 3, 0.0, R.data(), 3);

#endif

      // Now compute the SVD of R...
      char joba='G';
      char jobu = 'U', jobv = 'V';
      int mv = 0;
      constexpr f77int lwork = 100;
      f77int m = 3, lda = 3, ldv = 3, info;
      double work[lwork];
      f77int nn = 3;
      vecDouble S(3);
      vecDouble V(9);

      dgesvj_(&joba, &jobu, &jobv, &m, &nn, R.data(), &lda, S.data(), &mv, V.data(), &ldv, work, &lwork, &info);

      if (info > 0) {
        char op = 'E';
        double eps = dlamch_(&op);
        std::cerr << boost::format("Warning- SVD in kabschCore() failed to converge with info=%d and TOL=%e\n") %
          info % (eps * sqrt(3.0));
        for (uint i=0; i<6; ++i)
          std::cerr << boost::format("         work[%d] = %f\n") % (i+1) % work[i];
        std::cerr << "\tThis may happen periodically, causing the output to be 'mostly' orthogonal.\n";
        std::cerr << "\tThe residual is typically small and should not appreciably affect the resulting\n";
        std::cerr << "\tsuperposition.  If this warning appears frequently, then please notify the LOOS\n";
        std::cerr << "\tdevelopers at loos.maintainer@gmail.com\n";
        //  throw(NumericalError("SVD in alignment::kabschCore returned an error", info));
      } else if (info < 0)
        throw(NumericalError("SVD in alignment::kabschCore returned an error", info));


      double dR = R[0]*R[4]*R[8] + R[3]*R[7]*R[2] + R[6]*R[1]*R[5] -
        R[0]*R[7]*R[5] - R[3]*R[1]*R[8] - R[6]*R[4]*R[2];


      double dV = V[0]*V[4]*V[8] + V[3]*V[7]*V[2] + V[6]*V[1]*V[5] -
        V[0]*V[7]*V[5] - V[3]*V[1]*V[8] - V[6]*V[4]*V[2];


      if (dR * dV < 0.0) {
        S[2] = -S[2];

        R[6] = -R[6];
        R[7] = -R[7];
        R[8] = -R[8];
      }

      return SVDTupleVec(R, S, V);
    }


    GCoord centerAtOrigin(vecDouble& v) {
      GCoord c;

      for (uint i=0; i<v.size(); i += 3) {
        c.x() += v[i];
        c.y() += v[i+1];
        c.z() += v[i+2];
      }

      for (uint i=0; i<3; ++i)
        c[i] = 3*c[i]/v.size();

      for (uint i=0; i<v.size(); i += 3) {
        v[i] -= c.x();
        v[i+1] -= c.y();
        v[i+2] -= c.z();
      }

      return c;
    }



    // Return the RMSD only for a kabsch alignment between U and V assuming
    // both are centered
    double centeredRMSD(const vecDouble& U, const vecDouble& V) {

      int n = U.size();

      double ssu[3] = {0.0, 0.0, 0.0};
      double ssv[3] = {0.0, 0.0, 0.0};

      for (int j=0; j<n; j += 3) {
        for (uint i=0; i<3; ++i) {
          ssu[i] += U[j+i] * U[j+i];
          ssv[i] += V[j+i] * V[j+i];
        }
      }

      n /= 3;

      double E0 = ssu[0] + ssu[1] + ssu[2] + ssv[0] + ssv[1] + ssv[2];

      SVDTupleVec svd = kabschCore(U, V);

      vecDouble S(boost::get<1>(svd));
      double ss = S[0] + S[1] + S[2];
      double rmsd = std::sqrt(std::abs(E0-2.0*ss)/n);

      return(rmsd);
    }



    // Return the RMSD only for a kabsch alignment between U and V
    // Both will be centered first.
    double alignedRMSD(const vecDouble& U, const vecDouble& V) {

      int n = U.size();

      vecDouble cU(U);
      vecDouble cV(V);

      centerAtOrigin(cU);
      centerAtOrigin(cV);

      SVDTupleVec svd = kabschCore(cU, cV);

      double ssu[3] = {0.0, 0.0, 0.0};
      double ssv[3] = {0.0, 0.0, 0.0};

      for (int j=0; j<n; j += 3) {
        for (uint i=0; i<3; ++i) {
          ssu[i] += cU[j+i] * cU[j+i];
          ssv[i] += cV[j+i] * cV[j+i];
        }
      }

      n /= 3;

      double E0 = ssu[0] + ssu[1] + ssu[2] + ssv[0] + ssv[1] + ssv[2];

      vecDouble S(boost::get<1>(svd));
      double ss = S[0] + S[1] + S[2];
      double rmsd = std::sqrt(std::abs(E0-2.0*ss)/n);
      return(rmsd);
    }




    // Kabsch alignment between U and V, assuming both are centered.
    // Returns the tranformation matrix to align U onto V.
    GMatrix kabschCentered(const vecDouble& U, const vecDouble& V) {
      SVDTupleVec svd = kabschCore(U, V);

      vecDouble R(boost::get<0>(svd));
      vecDouble VV(boost::get<2>(svd));

      double M[9];

#if defined(__linux__) || defined(__CYGWIN__) || defined(__FreeBSD__)
      char ta = 'N';
      char tb = 'T';
      f77int three = 3;
      double one = 1.0;
      double zero = 0.0;

      dgemm_(&ta, &tb, &three, &three, &three, &one, R.data(), &three, VV.data(), &three, &zero, M, &three);

#else

      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 3, 3, 3, 1.0, R.data(), 3, VV.data(), 3, 0.0, M, 3);

#endif

      // Construct the new transformation matrix...  (W = M')
      GMatrix Z;
      for (uint i=0; i<3; i++)
        for (uint j=0; j<3; j++)
          Z(i,j) = M[i*3+j];

      return Z;

    }


    GMatrix kabsch(const vecDouble& U, const vecDouble& V) {

      vecDouble cU(U);
      vecDouble cV(V);

      GCoord U_center = centerAtOrigin(cU);
      GCoord V_center = centerAtOrigin(cV);
      GMatrix M = kabschCentered(cU, cV);

      XForm W;
      W.identity();
      W.translate(V_center);
      W.concat(M);
      W.translate(-U_center);

      return W.current();
    }


    void applyTransform(const GMatrix& M, vecDouble& v) {
      for (uint i=0; i<v.size(); i += 3) {
        GCoord c(v[i],v[i+1],v[i+2]);

        c = M * c;
        v[i] = c.x();
        v[i+1] = c.y();
        v[i+2] = c.z();
      }
    }


    vecDouble averageCoords(const vecMatrix& ensemble) {
      uint m = ensemble.size();
      uint n = ensemble[0].size();

      vecDouble avg(n, 0.0);
      for (uint j=0; j<m; ++j)
        for (uint i=0; i<n; ++i)
          avg[i] += ensemble[j][i];

      for (uint i=0; i<n; ++i)
        avg[i] /= m;

      return avg;
    }


    double rmsd(const vecDouble& u, const vecDouble& v) {
      double rms = 0.0;
      for (uint i=0; i<u.size(); i += 3) {
        double l = 0.0;
        for (uint j=0; j<3; ++j) {
          double d = u[i+j] - v[i+j];
          l += d*d;
        }
        rms += l;
      }
      rms = std::sqrt(3.0 * rms / u.size());

      return rms;
    }



  }


  // The following are the routines most should use...  They behave the same way as the old
  // ones from ensembles.cpp


  boost::tuple<std::vector<XForm>,greal,int> iterativeAlignment(alignment::vecMatrix& ensemble,
                                                                greal threshold, int maxiter) {
    using namespace alignment;

    int n = ensemble.size();
    std::vector<XForm> xforms(n);

    // Start by aligning against the first structure in the ensemble

    vecDouble target(ensemble[0]);
    centerAtOrigin(target);

    double rms;
    int iter = 0;

    do {
      for (int i = 0; i<n; i++) {
        XForm M(kabsch(ensemble[i], target));
        applyTransform(M.current(), ensemble[i]);
        xforms[i].premult(M.current());
      }

      vecDouble avg = averageCoords(ensemble);
      rms = rmsd(target, avg);
      target = avg;
      ++iter;
    } while (rms > threshold && iter <= maxiter );

    boost::tuple<std::vector<XForm>, greal, int> res(xforms, rms, iter);
    return(res);
  }


  boost::tuple<std::vector<XForm>,greal,int> iterativeAlignment(std::vector<AtomicGroup>& ensemble,
                                                                greal threshold, int maxiter) {
    using namespace alignment;

    int n = ensemble.size();
    std::vector<XForm> xforms(n);

    // Start by aligning against the first structure in the ensemble
    vecDouble target(ensemble[0].coordsAsVector());
    centerAtOrigin(target);

    double rms;
    int iter = 0;
    do {
      for (int i = 0; i<n; i++) {
        XForm M(kabsch(ensemble[i].coordsAsVector(), target));
        ensemble[i].applyTransform(M);
        xforms[i].premult(M.current());
      }

      AtomicGroup avg_structure = averageStructure(ensemble);
      vecDouble avg = avg_structure.coordsAsVector();
      rms = rmsd(target, avg);
      target = avg;
      ++iter;
    } while (rms > threshold && iter <= maxiter );

    boost::tuple<std::vector<XForm>, greal, int> res(xforms, rms, iter);
    return(res);
  }




  boost::tuple<std::vector<XForm>, greal, int> iterativeAlignment(const AtomicGroup& g,
                                                                  pTraj& traj,
                                                                  const std::vector<uint>& frame_indices,
                                                                  greal threshold, int maxiter) {

    using namespace alignment;

    // Must first prime the loop...
    AtomicGroup frame = g.copy();
    traj->readFrame(frame_indices[0]);
    traj->updateGroupCoords(frame);

    int nf = frame_indices.size();

    int iter = 0;
    greal rms;
    std::vector<XForm> xforms(nf);
    AtomicGroup avg = frame.copy();

    AtomicGroup target = frame.copy();
    target.centerAtOrigin();

    do {
      // Compute avg internally so we don't have to read traj twice...
      for (uint j=0; j<avg.size(); ++j)
        avg[j]->coords() = GCoord(0,0,0);

      for (int i=0; i<nf; ++i) {

        traj->readFrame(frame_indices[i]);
        traj->updateGroupCoords(frame);

        GMatrix M = frame.alignOnto(target);
        xforms[i].load(M);

        for (uint j=0; j<avg.size(); ++j)
          avg[j]->coords() += frame[j]->coords();
      }

      for (uint j=0; j<avg.size(); ++j)
        avg[j]->coords() /= nf;

      rms = avg.rmsd(target);
      target = avg.copy();
      ++iter;
    } while (rms > threshold && iter <= maxiter);

    boost::tuple<std::vector<XForm>, greal, int> res(xforms, rms, iter);
    return(res);
  }


  boost::tuple<std::vector<XForm>, greal, int> iterativeAlignment(const AtomicGroup& g,
                                                                  pTraj& traj,
                                                                  greal threshold, int maxiter) {

    std::vector<uint> framelist(traj->nframes());
    for (uint i=0; i<traj->nframes(); ++i)
      framelist[i] = i;

    return(iterativeAlignment(g, traj, framelist, threshold, maxiter));


  }

}
