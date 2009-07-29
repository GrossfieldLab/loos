#if !defined(XDRFILE_HPP)
#define XDRFILE_HPP

#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <exception>
#include <vector>

#include <loos_defs.hpp>
#include <rpc/xdr.h>


#if defined(__APPLE__)

typedef void* xdr_arg_t;

#else  /* defined(__APPLE__) */

typedef char* xdr_arg_t;

#endif  /* defined(__APPLE__) */




namespace loos {


  namespace internal {

    namespace xdr {

      extern "C" int overflow(xdr_arg_t p, xdr_arg_t data, int n);
      extern "C" int underflow(xdr_arg_t p, xdr_arg_t data, int n);

      int wrapper(XDR* x, int* i);
      int wrapper(XDR* x, float* f);
      int wrapper(XDR* x, double* d);
      int wrapper(XDR* x, unsigned int* u);
      int wrapper(XDR* x, long* l);
      int wrapper(XDR* x, unsigned long* u);
    }

    class XDRFile {
    public:

      XDRFile(std::iostream* is, std::ios_base::openmode mode = std::ios_base::in) : closed(false), stream(is) {
        xdrrec_create(&xdr, 0, 0, reinterpret_cast<xdr_arg_t>(stream), &xdr::underflow, &xdr::overflow);
        if (mode == std::ios_base::in) {
          xdr.x_op = XDR_DECODE;
          xdrrec_skiprecord(&xdr);
          writing = false;
        } else {
          xdr.x_op = XDR_ENCODE;
          writing = true;
        }
      }

      ~XDRFile() { close(); }

      void close(void) {
        if (closed)
          return;

        if (writing)
          xdrrec_endofrecord(&xdr, true);
        xdr_destroy(&xdr);
        closed = true;
      }

      std::iostream* operator()(void) { return(stream); }

      template<typename T> int read(std::vector<T>& data, const int n) {
        T datum;

        int i;
        for (i=0; i<n && xdr::wrapper(&xdr, &datum); ++i)
          data.push_back(datum);

        return(i);
      }

      template<typename T> int read(T& datum) {
        return(xdr::wrapper(&xdr, &datum));
      }

      template<typename T> int read(T* data, const int n) {
        int i;
        for (i=0; i<n && xdr::wrapper(&xdr, &(data[i])); ++i) ;
        return(i);
      }

      int read(char* data, const int n) {
        return(xdr_opaque(&xdr, data, n));
      }


      template<typename T> int write(std::vector<T>& data) {
        int n=0;
        for (typename std::vector<T>::iterator i(data.begin());
             i != data.end() && xdr::wrapper(&xdr, &(*i));
             ++i, ++n)
          ;
        return(n);
      }

      template<typename T> int write(T& datum) {
        return(xdr::wrapper(&xdr, &datum));
      }

      template<typename T> int write(T* data, const int n) {
        int i;
        for (i=0; i<n && xdr::wrapper(&xdr, &(data[i])); ++i) ;
        return(i);
      }
          
      int write(char* data, const int n) {
        return(xdr_opaque(&xdr, data, n));
      }

    private:
      bool writing, closed;
      std::iostream* stream;
      XDR xdr;
    };

  } /* internal */
} /* loos */



#endif
