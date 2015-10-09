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


        }

#if !defined(SWIG)
        boost::tuple<std::vector<XForm>,greal,int> iterativeAlignment(alignment::vecMatrix& ensemble,
                                                                      greal threshold=1e-6,
                                                                      int maxiter=1000);
        
        boost::tuple<std::vector<XForm>,greal,int> iterativeAlignment(std::vector<AtomicGroup>& ensemble,
                                                                      greal threshold=1e-6,
                                                                      int maxiter=1000);

        //! Compute an iterative superposition by reading in frames from the Trajectory.
        /**
         * The iterativeAlignment() functions that take a trajectory as an argument do
         * NOT cache frames of the trajectory internally.  This means that the trajectory
         * will be read as many times as is necessary for the alignment to
         * converge.  In practice, the OS-specific caching will likely result
         * in decent performance.  If speed is essential, then consider
         * using the iterativeAlignment() version that takes a
         * \p std::vector<AtomicGroup>& as argument instead.
         */
        boost::tuple<std::vector<XForm>,greal,int> iterativeAlignment(const AtomicGroup& model,
                                                                      pTraj& traj,
                                                                      const std::vector<uint>& frame_indices,
                                                                      greal threshold=1e-6,
                                                                      int maxiter=1000);


        boost::tuple<std::vector<XForm>,greal,int> iterativeAlignment(const AtomicGroup& model,
                                                                      pTraj& traj,
                                                                      greal threshold=1e-6,
                                                                      int maxiter=1000);

        
#endif // !defined(SWIG)
}



#endif

    

