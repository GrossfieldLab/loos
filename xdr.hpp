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
    public:
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

      uint read(char** p) {
        uint n;

        if (!read(n))
          return(0);
        char* s = new char[n+1];
        s[n] = '\0';
        uint i = read(s, n);
        *p = s;
        return(i);
      }

      uint read(std::string& s) {
        char* p;
        uint i = read(&p);
        if (!i)
          return(0);

        s = std::string(p);
        return(i);
      }


      // -----------------------------------------------------

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

      template<typename T> uint write(const T* ary, const uint n) {
        uint i;
        for (i=0; i<n && write(&(ary[i])); ++i) ;
        return(i);
      }

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
