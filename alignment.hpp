#if !defined(ALIGNMENT_HPP)
#define ALIGNMENT_HPP

#include <vector>
#include <boost/tuple/tuple.hpp>

#include <loos_defs.hpp>
#include <MatrixImpl.hpp>
#include <MatrixOps.hpp>

#include <XForm.hpp>


namespace loos {
    
    namespace alignment {
        
        typedef std::vector<double>    vecDouble;
        typedef std::vector<vecDouble> vecMatrix;
        typedef boost::tuple<vecDouble, vecDouble, vecDouble>   SVDTupleVec;
        
        
        SVDTupleVec kabschCore(const vecDouble& u, const vecDouble& v);
        GCoord centerAtOrigin(vecDouble& v);
        double alignedRMSD(const vecDouble& U, const vecDouble& V);
        double centeredRMSD(const vecDouble& U, const vecDouble& V);
        GMatrix kabsch(const vecDouble& U, const vecDouble& V);
        void applyTransform(const GMatrix& M, vecDouble& v);
        vecDouble averageCoords(const vecMatrix& ensemble);
        double rmsd(const vecDouble& u, const vecDouble& v);
        
        boost::tuple<std::vector<XForm>,greal,int> iterativeAlignment(vecMatrix& ensemble,
                                                                      greal threshold=1e-6,
                                                                      int maxiter=1000);
        
        boost::tuple<std::vector<XForm>,greal,int> iterativeAlignment(std::vector<AtomicGroup>& ensemble,
                                                                      greal threshold=1e-6,
                                                                      int maxiter=1000);
        
        boost::tuple<std::vector<XForm>,greal,int> iterativeAlignment(const AtomicGroup& model,
                                                                      pTraj& traj,
                                                                      const std::vector<uint>& frame_indices,
                                                                      greal threshold=1e-6,
                                                                      int maxiter=1000);

  }

}



#endif

    

