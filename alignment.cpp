#include <loos_defs.hpp>
#include <alignment.hpp>

namespace loos {



  namespace alignment {

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
      f77int m = 3, lda = 3, ldv = 3, lwork=100, info;
      double work[lwork];
      f77int nn = 3;
      vecDouble S(3);
      vecDouble V(9);
  
      dgesvj_(&joba, &jobu, &jobv, &m, &nn, R.data(), &lda, S.data(), &mv, V.data(), &ldv, work, &lwork, &info);
    
      if (info != 0)
	throw(NumericalError("SVD in AtomicGroup::superposition returned an error", info));


  
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



    double alignedRMSD(const vecDouble& U, const vecDouble& V) {

      int n = U.size()/3;

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
      double rmsd = sqrt(abs(E0-2.0*ss)/n);
      return(rmsd);
    }


  
  

    XForm kabsch(const vecDouble& U, const vecDouble& V) {
    
      vecDouble cU(U);
      vecDouble cV(V);

      GCoord U_center = centerAtOrigin(cU);
      GCoord V_center = centerAtOrigin(cV);


      SVDTupleVec svd = kabschCore(cU, cV);

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

      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, R.data(), 3, VV.data(), 3, 0.0, M, 3);

#endif

      // Construct the new transformation matrix...  (W = M')
      GMatrix Z;
      for (uint i=0; i<3; i++)
	for (uint j=0; j<3; j++)
	  Z(i,j) = M[i*3+j];

      XForm W;
      W.identity();
      W.translate(V_center);
      W.concat(Z);
      W.translate(-U_center);

      return W;
    }


    // hack for now...
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
      rms = sqrt(3.0 * rms / u.size());

      return rms;
    }


    boost::tuple<std::vector<XForm>,greal,int> iterativeAlignment(vecMatrix& ensemble,
								  greal threshold, int maxiter) {
      int n = ensemble.size();
      std::vector<XForm> xforms(n);

      // Start by aligning against the first structure in the ensemble
      vecDouble target(ensemble[0]);
      centerAtOrigin(target);
      for (uint i=1; i<n; ++i) {
	XForm M = kabsch(ensemble[i], target);
	applyTransform(M.current(), ensemble[i]);
	xforms[i].premult(M.current());
      }
      target = averageCoords(ensemble);

      double rms;
      uint iter = 0;
    
      do {
	for (int i = 0; i<n; i++) {
	  XForm M = kabsch(ensemble[i], target);
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



  }



  
}
