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


#if !defined(LOOS_TRAJWRITER_HPP)
#define LOOS_TRAJWRITER_HPP

#include <iostream>
#include <string>
#include <stdexcept>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include <loos_defs.hpp>
#include <AtomicGroup.hpp>


namespace loos {


  //! Base class for writing trajectories
  /**
   * This is the interface for creating and writing to trajectories.
   * The interface is kept simple on purpose to make it easy to work
   * with multiple formats.  This means that trajectory and frame
   * metadata may be set to default values.  If you need more control
   * over how the trajectory is written, then use the derived classes
   * explicitly.
   */

  class TrajectoryWriter {
  public:

    //! Exception while writing
    struct WriteError : public std::exception {
      WriteError() : text("Error while writing trajectory") {}
      WriteError(const char* message) : text(message) {}
      
      virtual const char* what() const throw() { return(text); }
      
      const char* text;
    };
    

    //! Write a trajectory to a file, optionally appending
    TrajectoryWriter(const std::string& fname, const bool append = false)
      : appending_(false) {
      struct stat statbuf;

      if (append && !stat(fname.c_str(), &statbuf))
	openStream(fname, true);
      else
	openStream(fname);
    }


    //! Write a trajectory to a stream
    /**
     * Note that this constructor assumes that the stream is correctly
     * prepped by the caller...if you need to seekp() to the end of the
     * stream for appending, then this must be done before instantiating
     * the TrajectoryWriter object.
     */
    TrajectoryWriter(std::iostream* s, const bool append = false)
      : appending_(append), delete_(false) {}


    virtual ~TrajectoryWriter() {
      if (delete_)
	delete stream_;
    }


    //! Set comments in metadata (not all formats support)
    virtual void setComments(const std::vector<std::string>& comments) { }

    //! Set comment in metadata (not all formats support)
    virtual void setComments(const std::string& s) {
      std::vector<std::string> c(1, s);
      setComments(c);
    }

    //! Wirte a single frame
    virtual void writeFrame(const AtomicGroup& model) =0;

    //! Write a single frame specifying the step and timepoint
    /**
     * Not all formats support this.  By default, it will drop the
     * extra data and call writeFrame()
     */
    virtual void writeFrame(const AtomicGroup& model, const uint step, const double time) {
      writeFrame(model);
    }

    //! Can format write step on a per-frame basis?
    virtual bool hasFrameStep() const { return(false); }

    //! Can format write time on a per-frame basis?
    virtual bool hasFrameTime() const { return(false); }

    //! Does format support comments in metadata?
    virtual bool hasComments() const { return(false); }

    //! Total frames in output file
    /**
     * For files being appended too, this includes the frames already
     * written...
     */
    virtual uint framesWritten() const =0;

    //! Returns true if appending to an existing trajectory
    bool isAppending() const { return(appending_); }

  protected:
    std::iostream* stream_;
    bool appending_;
    bool delete_;

  private:


    // Handle opening up a stream to a file...  If it exists and we
    // are asked to append, seek to the end of the file.
    
    void openStream(const std::string& fname, const bool append = false) {

      

      std::ios_base::openmode mode = std::ios_base::out | std::ios_base::binary;
      if (append)
        mode |= std::ios_base::in;
      else
        mode |= std::ios_base::trunc;

      stream_ = new std::fstream(fname.c_str(), mode);
      if (append) {
	stream_->seekp(0, std::ios_base::end);
	// Check to see if file is empty...
	if (stream_->tellp() == 0)
	  appending_ = false;
	else
	  appending_ = true;
      }

      if (!stream_->good())
        throw(std::runtime_error("Error while opening output trajectory file"));


      delete_ = true;    // Delete the stream pointer when dtor called
    }
    

  };



}




#endif
