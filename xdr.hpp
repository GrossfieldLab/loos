/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo, Alan Grossfield
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


#if !defined(XDR_HPP)
#define XDR_HPP

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include <loos_defs.hpp>
#include <utils.hpp>


namespace loos {

  namespace internal {

    //! This class provides some facility for handling XDR data
    /**
     * The read and write functions use templates to read the
     * appropriate raw data.  Beware of unexpected type
     * conversions...  All functions also return a 0 for error, or a 1
     * for success (or the number of elements actually read/written).
     */
    class XDR {
    public:
      //! Type (and hence size) of the external block
      typedef unsigned int           block_type;

    public:

      //! Constructor determines need to convert data at instantiation
      XDR(std::iostream* s) : stream(s), need_to_swab(false) {
        int test = 0x1234;
        if (*(reinterpret_cast<char*>(&test)) == 0x34) {
          need_to_swab = true;
        }
      }

      //! Returns the stored iostream pointer
      std::iostream* get(void) { return(stream); }

      //! Read a single datum
      template<typename T> uint read(T* p) {

        if (sizeof(T) > sizeof(block_type))
          throw(std::logic_error("Attempting to read a POD that is too large"));

        T result;
        stream->read(reinterpret_cast<char*>(&result), sizeof(block_type));
        if (sizeof(T) > 1 && need_to_swab)
          result = swab(result);

        *p = result;
        return(!stream->fail());
      }

      template<typename T> uint read(T& t) { return(read(&t)); }

      //! Read an n-array of data
      template<typename T> uint read(T* ary, const uint n) {
        uint i;
        for (i=0; i<n && read(ary+i); ++i) ;
        return(i);
      }


      //! Read in an opaque array of n-bytes (same as xdr_opaque)
      uint read(char* p, uint n) {
        uint rndup;
        static char buf[sizeof(block_type)];

        if (n == 0)
          return(1);

        rndup = n % sizeof(block_type);
        if (rndup > 0)
          rndup = sizeof(block_type) - rndup;

        stream->read(p, n);
        if (stream->fail())
          return(0);
        if (rndup)
          stream->read(buf, rndup);

        return(n);
      }

      //! Same as xdr_string.

      uint read(boost::shared_ptr<char>& p) {
        uint n;

        if (!read(n))
          return(0);
        char* s = new char[n+1];
        uint i = read(s, n);
        s[n] = '\0';
        p = boost::shared_ptr<char>(s);
        return(i);
      }

      uint read(std::string& s) {
        boost::shared_ptr<char> p;
        int i = read(p);
        if (!i)
          return(0);

        s = std::string(p.get());

        return(i);
      }


      // -----------------------------------------------------

      //! Writes a single datum
      template<typename T> uint write(const T* p) {

        if (sizeof(T) > sizeof(block_type))
          throw(std::logic_error("Attempting to write a POD that is too large"));

        uint u;
        T* up = reinterpret_cast<T*>(&u);
        *up = *p;

        if (sizeof(T) > 1 && need_to_swab)
          u = swab(u);

    
        stream->write(reinterpret_cast<char*>(&u), sizeof(block_type));

        return(!stream->fail());
      }

      template<typename T> uint write(const T& t) { return(write(&t)); }

      //! Writes an n-array of data
      template<typename T> uint write(const T* ary, const uint n) {
        uint i;
        for (i=0; i<n && write(&(ary[i])); ++i) ;
        return(i);
      }

      //! Writes an opaque array of n-bytes
      uint write(const char* p, const uint n) {
        uint rndup;
        static char buf[sizeof(block_type)];
        static bool init(false);

        if (!init)
          for (uint i=0; i<sizeof(block_type); ++i)
            buf[i] = '\0';

        rndup = n % sizeof(block_type);
        if (rndup > 0)
          rndup = sizeof(block_type) - rndup;

        stream->write(p, n);
        if (!stream->fail())
          stream->write(buf, rndup);

        return(stream->fail() ? 0 : n);
      }

      //! Writes a C-string (ie xdr_string)
      uint write(const char* p) {
        uint n = strlen(p);
        write(n);
        return(write(p, n));
      }

      uint write(const std::string& s) { return(write(s.c_str())); }

    private:
      std::iostream* stream;
      bool need_to_swab;

    };




  } /* internal */
} /* loos */



#endif
