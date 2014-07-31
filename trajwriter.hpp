#if !defined(LOOS_TRAJWRITER_HPP)
#define LOOS_TRAJWRITER_HPP

#include <iostream>
#include <string>
#include <stdexcept>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>


namespace loos {

  class TrajectoryWriter {
  public:

    struct WriteError : public std::exception {
      WriteError() : text("Error while writing trajectory") {}
      WriteError(const char* message) : text(message) {}
      
      virtual const char* what() const throw() { return(text); }
      
      const char* text;
    };
    
    virtual void writeFrame(const AtomicGroup& model) =0;
  };



}




#endif
