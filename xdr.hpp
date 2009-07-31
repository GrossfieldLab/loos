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


    class XDR {
      typedef unsigned int           block_type;
    public:
      XDR(std::iostream* s) : stream(s), need_to_swab(false) {
        int test = 0x1234;
        if (*(reinterpret_cast<char*>(&test)) == 0x34) {
          need_to_swab = true;
        }
      }


      std::iostream* get(void) { return(stream); }


      template<typename T> uint read(T* p) {

        if (sizeof(T) > sizeof(block_type))
          throw(std::logic_error("Attempting to read a POD that is too large"));

        uint data;
        stream->read(reinterpret_cast<char*>(&data), sizeof(block_type));

        T* pdata = reinterpret_cast<T*>(&data);
        T result(*pdata);
        if (sizeof(T) > 1 && need_to_swab)
          result = swab(result);

        *p = result;
        return(!stream->fail());
      }

      template<typename T> uint read(T& t) { return(read(&t)); }

      template<typename T> uint read(T* ary, const uint n) {
        uint i;
        for (i=0; i<n && read(&(ary[i])); ++i) ;
        return(i);
      }


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


      // -----------------------------------------------------

      template<typename T> uint write(T* p) {

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

      template<typename T> uint write(T& t) { return(write(&t)); }

      template<typename T> uint write(T* ary, const uint n) {
        uint i;
        for (i=0; i<n && write(&(ary[i])); ++i) ;
        return(i);
      }

      uint write(char* p, const uint n) {
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

    private:
      std::iostream* stream;
      bool need_to_swab;

    };




  } /* internal */
} /* loos */



#endif
