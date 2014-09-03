/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2014, Tod D. Romo, Alan Grossfield
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


#if !defined(LOOS_XTCWRITER_HPP)
#define LOOS_XTCWRITER_HPP

#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>
#include <xdr.hpp>
#include <trajwriter.hpp>

namespace loos {

  //! Class for writing Gromacs XTC trajectories
  /**
   * This code borrows heavily from the xdrfile-1.1b library provided
   * by Gromacs.  By default, the writer will assume that the frames
   * are evenly spaced and will use the dt and steps_per_frame
   * variables to determine the step and timepoint for each frame.
   * For non-uniform frame intervals, explicitly pass a step and time
   * to writeFrame().  Note that this will NOT modify the internal
   * counters, so you should use on form of writeFrame() or the other
   * and not mix them.  If you must, use currentStep() to update the
   * internal step counter (and possibly timePerStep()).
   */


  class XTCWriter : public TrajectoryWriter {
    static const int magicints[];
    static const int firstidx;
    
    static const int lastidx;
    
    static const int DIM;

  public:

    static const float default_precision = 1e3;


    struct InternalOverflow : public WriteError { virtual const char* what() const throw() { return("Internal overflow compressing coordinates"); } };
    

    //! Class factory function
    static pTrajectoryWriter create(const std::string& s, const bool append = false) {
      return(pTrajectoryWriter(new XTCWriter(s, append)));
    }


    XTCWriter(const std::string fname, const bool append = false) :
      TrajectoryWriter(fname, append),
      buf1size(0), buf2size(0),
      buf1(0), buf2(0),
      natoms_(0),
      dt_(1.0),
      step_(0),
      steps_per_frame_(1),
      current_(0),
      crds_size_(0),
      crds_(0),
      precision_(default_precision)
    {
      xdr.setStream(stream_);
      if (appending_)
	prepareToAppend();
    }


    XTCWriter(const std::string fname, const double dt, const uint steps_per_frame, const bool append = false) :
      TrajectoryWriter(fname, append),
      buf1size(0), buf2size(0),
      buf1(0), buf2(0),
      natoms_(0),
      dt_(dt),
      step_(0),
      steps_per_frame_(steps_per_frame),
      current_(0),
      crds_size_(0),
      crds_(0),
      precision_(default_precision)
    {
      xdr.setStream(stream_);
      if (appending_)
	prepareToAppend();
    }
    



    ~XTCWriter() {
      delete[] buf1;
      delete[] buf2;
      delete[] crds_;
    }


    //! Get the time per step
    double timePerStep() const { return(dt_); }

    //! Set the time per step
    void timePerStep(const double dt) { dt_ = dt; }

    //! How many steps per frame written
    uint stepsPerFrame() const { return(steps_per_frame_); }

    //! Set how many steps pass per frame written
    void stepsPerFrame(const uint s) { steps_per_frame_ = s; }

    //! What the current output step is
    uint currentStep() const { return(step_); }

    //! Sets the current output step
    void currentStep(const uint s) { step_ = s; }


    //! Write a frame to the trajectory
    void writeFrame(const AtomicGroup& model);

    //! Write a frame to the trajectory with explicit step and time metadata
    void writeFrame(const AtomicGroup& model, const uint step, const double time);

    uint framesWritten() const { return(current_); }

  private:
    int sizeofint(const int size) const;
    int sizeofints(const int num_of_bits, const unsigned int sizes[]) const;
    void encodebits(int* buf, int num_of_bits, const int num) const;
    void encodeints(int* buf, const int num_of_ints, const int num_of_bits,
		    const unsigned int* sizes, const unsigned int* nums) const;
    void writeCompressedCoordsFloat(float* ptr, int size, float precision);
       
    void allocateBuffers(const size_t size);

    void writeHeader(const int natoms, const int step, const float time);
    void writeBox(const GCoord& box);

    void prepareToAppend();
    
  private:
    uint buf1size, buf2size;
    int* buf1;
    int* buf2;
    uint natoms_;
    double dt_;
    uint step_;
    uint steps_per_frame_;
    uint current_;
    uint crds_size_;
    float* crds_;
    float precision_;

    internal::XDRWriter xdr;
  };


};



#endif
