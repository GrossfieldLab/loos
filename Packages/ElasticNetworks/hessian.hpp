/*
  hessian matrix construction code

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


/** \addtogroup ENM
 *@{
 */



#if !defined(HESSIAN_HPP)
#define HESSIAN_HPP


#include <loos.hpp>

#include "spring_functions.hpp"

namespace ENM {

  //! This class creates superblocks in a hessian
  /**
   *This is the base class for defining elements (superblocks) in the
   * Hessian.  Note that this class does NOT use NVI, so derived classes should
   * verify that they have valid args for the block() function.
   */
  class SuperBlock {
  public:
    SuperBlock() : springs(0) { }
  
    //! Constructor taking a spring function and a list of nodes
    /**
     * Arguments:
     * \arg \c func SpringFunction to use to determine spring constants between nodes
     * \arg nodelist List of nodes in the model
     *
     * For example, to compute a hessian using the HCA method,
     \code
     SpringFunction *spring = new HCA;
     SuperBlock *blocker = new SuperBlock(spring, model);
     \endcode
    */
    SuperBlock(SpringFunction* func, const loos::AtomicGroup& nodelist) : springs(func), nodes(nodelist) { }
    SuperBlock(const SuperBlock& b) : springs(b.springs), nodes(b.nodes) { }
    virtual ~SuperBlock() { }

    uint size() const { return(static_cast<uint>(nodes.size())); }

    // ------------------------------------------------------
    //! Forwards to the contained SpringFunction...
    virtual SpringFunction::Params setParams(const SpringFunction::Params& v) {
      return(springs->setParams(v));
    }

    //! Forwards to the contained SpringFunction...
    virtual bool validParams() const { return(springs->validParams()); }

    //! Forwards to the contained SpringFunction...
    virtual uint paramSize() const { return(springs->paramSize()); }
    // ------------------------------------------------------

    //! Returns a 3x3 matrix representing a superblock in the Hessian for the two nodes
    virtual loos::DoubleMatrix block(const uint j, const uint i) {
      return(blockImpl(j, i, springs));
    }


  protected:

    //! Implementation of the superblock calculation
    /**
     *This is the actual implementation of the SuperBlock calculation.
     *In most cases, derived clases will probably want to use this but
     * with alternative spring functions...
     */
    loos::DoubleMatrix blockImpl(const uint j, const uint i, SpringFunction* fptr) {
      if (i >= size() || j >= size())
        throw(std::runtime_error("Invalid index in Hessian SuperBlock"));

      if (fptr == 0)
        throw(std::runtime_error("No spring function defined for hessian!"));

      loos::GCoord u = nodes[j]->coords();
      loos::GCoord v = nodes[i]->coords();
      loos::GCoord d = v - u;
    
      loos::DoubleMatrix K = fptr->constant(u, v, d);
      loos::DoubleMatrix B(3, 3);
      for (uint y=0; y<3; ++y)
        for (uint x=0; x<3; ++x)
          B(x, y) = d[x]*d[y] * K(x,y);

      return(B);
    }


    SpringFunction* springs;
    loos::AtomicGroup nodes;
  };



  //! SuperBlock decorator base class
  /**
   *The following is a decorator for a SuperBlock.
   *It both inherits from (so it can be used in place of a SuperBlock)
   *and contains a SuperBlock.  This allows additional behavior to be
   *layed on top of the SuperBlock.
   */
  class SuperBlockDecorator : public SuperBlock {
  public:

    //! Constructor that takes a SuperBlock to decorate
    SuperBlockDecorator(SuperBlock* b) : SuperBlock(*b), decorated(b) { }

  protected:
    SuperBlock *decorated;
  };



  //! Decorator for switching spring functions based on a matrix of flags
  /**
   *The following is a decorator for SuperBlock that implements an
   *alternative set of spring constants for nodes that are "bound"
   *together.  The constructor takes a SuperBlock to decorate, along
   *with a pointer to the alternative SpringFunction and a matrix of
   *ints representing the connectivity (i.e. 1 if two nodes are
   *connected, 0 otherwise).
   *
   *A few notes about using decorators...  The idea behind a decorator
   *is that you add layers (or decorate) to a class by combining
   *multiple decorators.  For example, suppose you have two different
   *kinds of connectivity you want to represent in a Hessian.  You
   *would set-up your SuperBlock like,
   \code
   SuperBlock* unbound = new SuperBlock(unbound_spring, nodes);
   BoundSuperBlock* backboned = new BoundSuperBlock(unbound, backbone_springs, backbone_bonds);
   BoundSuperBlock* side_chained = new BoundSuperBlock(backboned, side_chain_springs, side_chain_bonds);
   \endcode
   *You now always work with the last decorated object,
   *i.e. side_chained.  When side_chained->block() is called, it first
   *checks to see if the nodes represent a side-chain bond.  If so,
   *that spring function is used.  If not, then it passes control to
   *the object it decorates, i.e. backboned.  Backboned now checks to
   *see if the nodes represent a backbone bond.  If so, it uses that
   *spring function.  If not, then control is passed to the inner
   *unbound SuperBlock which uses its spring function.
   *
   *This method has two important caveats.  First, the calculation is
   *now order-dependent.  If, for some reason, you have nodes that are
   *listed as both side-chains and backbones (for a contrived example),
   *then the one used will depend on the order in which the SuperBlock
   *was decorated.  The second caveat is that we are using real, raw
   *pointers here, so be careful about cleaning up to avoid memory
   *leaks and also keep in mind that the intermediate pointers
   *(i.e. backboned) are contained within the higher-level decorators.
   *So, do NOT delete any of the intermediate steps until you are sure
   *you are done with everything.
   */
  class BoundSuperBlock : public SuperBlockDecorator {
  public:
    BoundSuperBlock(SuperBlock* b, SpringFunction* bs, loos::Math::Matrix<int>& cm) :
      SuperBlockDecorator(b),
      bound_spring(bs),
      connectivity(cm)
    {
      if (connectivity.rows() != connectivity.cols() && connectivity.cols() != size())
        throw(std::runtime_error("Connectivity matrix and Nodelist have differing sizes"));
    }


    // Block now checks to see if nodes i and j are connected and, if
    // so, uses our alternative spring function.
    loos::DoubleMatrix block(const uint j, const uint i) {
      if (connectivity(j, i))
        return(blockImpl(j, i, bound_spring));
      else
        return(decorated->block(j, i));
    }

    //! Assign parameters and propagate to the decorated superblock
    SpringFunction::Params setParams(const SpringFunction::Params& v) {
      SpringFunction::Params u = bound_spring->setParams(v);
      if (! u.empty())
        u = decorated->setParams(u);
      return(u);
    }

    //! Both the alternate and all decorated parameters are valid
    bool validParams() const { return(bound_spring->validParams() && decorated->validParams()); }

    //! Returns the aggregate parameter size
    uint paramSize() const { return(bound_spring->paramSize() + decorated->paramSize()); }

  private:
    SpringFunction* bound_spring;
    loos::Math::Matrix<int> connectivity;
  };

};

#endif



/** @} */
