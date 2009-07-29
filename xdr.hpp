#if !defined(XDRFILE_HPP)
#define XDRFILE_HPP

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>


#include <loos_defs.hpp>
#include <utils.hpp>


namespace loos {

  namespace internal {


    class XDR {
    public:
      XDR(iostream* s) : stream(s), need_to_swab(false) {
        int test = 0x1234;
        if (*(reinterpret_cast<char*>(&test)) == 0x34) {
          need_to_swab = true;
        }
      }


      template<typename T> int read(T* p) {

        if (sizeof(T) > sizeof(uint))
          throw(logic_error("Attempting to read a POD that is too large"));

        uint data;
        stream->read(reinterpret_cast<char*>(&data), sizeof(uint));

        T* pdata = reinterpret_cast<T*>(&data);
        T result(*pdata);
        if (sizeof(T) > 1 && need_to_swab)
          result = swab(result);

        *p = result;
        return(!stream->fail());
      }

      template<typename T> int read(T& t) { return(read(&t)); }

      template<typename T> int read(T* ary, const int n) {
        int i;
        for (i=0; i<n && read(&(ary[i])); ++i) ;
        return(i);
      }


      int read(char* p, uint n) {
        uint rndup;
        static char buf[sizeof(uint)];

        rndup = n % sizeof(uint);
        if (rndup > 0)
          rndup = sizeof(uint) - rndup;

        stream->read(p, n);
        if (!stream->fail())
          stream->read(buf, rndup);

        return(stream->gcount());
      }


      // -----------------------------------------------------

      template<typename T> int write(T* p) {

        if (sizeof(T) > sizeof(uint))
          throw(logic_error("Attempting to write a POD that is too large"));

        uint u;
        T* up = reinterpret_cast<T*>(&u);
        *up = *p;

        if (sizeof(T) > 1 && need_to_swab)
          u = swab(u);

    
        stream->write(reinterpret_cast<char*>(&u), sizeof(uint));

        return(!stream->fail());
      }

      template<typename T> int write(T& t) { return(write(&t)); }

      template<typename T> int write(T* ary, const int n) {
        int i;
        for (i=0; i<n && write(&(ary[i])); ++i) ;
        return(i);
      }

      int write(char* p, uint n) {
        uint rndup;
        static char buf[sizeof(uint)];
        static bool init(false);

        if (!init)
          for (int i=0; i<sizeof(uint); ++i)
            buf[i] = '\0';

        rndup = n % sizeof(uint);
        if (rndup > 0)
          rndup = sizeof(uint) - rndup;

        stream->write(p, n);
        if (!stream->fail())
          stream->write(buf, rndup);

        return(stream->fail() ? 0 : n);
      }

    private:
      iostream* stream;
      bool need_to_swab;

    };




  } /* internal */
} /* loos */



#endif
