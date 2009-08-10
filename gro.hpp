#if !defined(GRO_HPP)
#define GRO_HPP


#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>

namespace loos {

  class Gromacs : public AtomicGroup {
  public:

    Gromacs() { }
    
    explicit Gromacs(const char* fname) {
      std::ifstream ifs(fname);
      if (!ifs)
        throw(std::runtime_error("Cannot open Gromacs file " + std::string(fname)));
      read(ifs);
    }

    explicit Gromacs(const std::string& fname) {
      std::ifstream ifs(fname.c_str());
      if (!ifs)
        throw(std::runtime_error("Cannot open Gromacs file " + fname));
      read(ifs);
    }

    explicit Gromacs(std::istream& ifs) { read(ifs); }


    std::string title(void) const { return(title_); }

  private:
    std::string title_;

  private:
    
    void read(std::istream& ifs);

  };



}


#endif
