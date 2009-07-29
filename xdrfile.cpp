#include <xdrfile.hpp>



namespace loos {


  namespace internal {

    namespace xdr {

      extern "C" int overflow(xdr_arg_t p, xdr_arg_t data, int n) {
        std::iostream* of(reinterpret_cast<std::iostream*>(p));
        of->write(static_cast<char*>(data), n);
        return( (of->fail() || of->bad()) ? -1 : n);
      }
      
      
      extern "C" int underflow(xdr_arg_t p, xdr_arg_t data, int n) {
        std::iostream* s(reinterpret_cast<std::iostream*>(p));
        
        if (s->eof())
          return(-1);
        s->read(static_cast<char*>(data), n);
        return( s->fail() ? s->gcount() : -1);
      }

      int wrapper(XDR* x, int* i) { return(xdr_int(x, i)); }
      int wrapper(XDR* x, float* f) { return(xdr_float(x, f)); }
      int wrapper(XDR* x, double* d) { return(xdr_double(x, d)); }
      int wrapper(XDR* x, unsigned int* u) { return(xdr_u_int(x, u)); }
      int wrapper(XDR* x, long* l) { return(xdr_long(x, l)); }
      int wrapper(XDR* x, unsigned long* u) { return(xdr_u_long(x, u)); }
      
    } /* internal */
  } /* loos */
}

