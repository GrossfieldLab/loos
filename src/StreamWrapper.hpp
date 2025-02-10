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

#if !defined(LOOS_STREAMWRAPPER_HPP)
#define LOOS_STREAMWRAPPER_HPP

#include <istream>
#include <fstream>
#include <string>
#include <stdexcept>

#include <boost/utility.hpp>
#include <exceptions.hpp>

namespace loos
{

  //! Simple wrapper class for caching stream pointers
  /** This class was written primarily for use with the DCD classes
   *  where we want to have a cached stream that we may read from (or
   *  write to) at various times in the future.  Access to the
   *  underlying fstream pointer is through the operator() functor.
   *
   *  The basic idea here is that you pass the class either a string or
   *  a char array and that will be opened into a new stream (reading by
   *  default).  The fstream pointer will be stored and when the wrapper
   *  object is destroyed, the stream is released & deleted.  If you
   *  pass the wrapper an fstream, however, the internal pointer is
   *  initialized to point to that stream and when the wrapper object is
   *  destroyed, the stream is left alone.
   */
  class StreamWrapper : public boost::noncopyable
  {
  public:
    StreamWrapper() : new_stream(false), stream(0) {}

    //! Sets the internal stream pointer to fs
    explicit StreamWrapper(std::istream &fs) : new_stream(false), stream(&fs) {}

    //! Opens a new stream with file named 's'
    StreamWrapper(const std::string &s,
                  const std::ios_base::openmode mode = std::ios_base::in | std::ios_base::binary)
        : new_stream(true)
    {
      stream = new std::fstream(s.c_str(), mode);
      if (!stream->good())
        throw(FileOpenError(s));
    }

    //! Sets the internal stream to point to a newly opened filed...
    void setStream(const std::string &s,
                   const std::ios_base::openmode mode = std::ios_base::in | std::ios_base::binary)
    {
      if (new_stream)
        delete stream;

      new_stream = true;
      stream = new std::fstream(s.c_str(), mode);
      if (!stream->good())
        throw(FileOpenError(s));
    }

    //! Sets the internal stream to the passed fstream.
    void setStream(std::istream &fs)
    {
      if (new_stream)
        delete stream;
      new_stream = false;
      stream = &fs;
    }

    //! Returns the internal istream pointer
    std::istream *operator()(void)
    {
      if (stream == 0)
        throw(LOOSError("Attempting to access an unset stream"));
      return (stream);
    }

    //! Returns true if the internal stream pointer is unset
    bool isUnset(void) const { return (stream == 0); }

    ~StreamWrapper()
    {
      if (new_stream)
        delete stream;
    }

  private:
    bool new_stream;
    std::istream *stream;
  };

}

#endif
