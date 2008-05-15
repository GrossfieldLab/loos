/*
  ParserSelector.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

*/



#if !defined(PARSERSELECTOR_HPP)
#define PARSERSELECTOR_HPP

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <tr1/memory>
#include <boost/tuple/tuple.hpp>

using namespace std;
using namespace tr1;

#include <loos.hpp>
#include <AtomicGroup.hpp>
#include <Kernel.hpp>
#include <Tokenizer.hpp>
#include <Parser.hpp>



class ParserSelector : public AtomSelector {
private:
  loos::Kernel kernel;

public:
  ParserSelector(const string s) {
    bool b = loos::parser::parse(s, kernel);
    if (!b)
      throw(runtime_error("Unable to parse " + s));
  }

  bool operator()(const pAtom& pa) {
    kernel.execute(pa);
    if (kernel.stack().size() != 1)
      throw(runtime_error("Execution error - unexpected values on stack"));

    loos::Value results = kernel.stack().pop();
    if (results.type != loos::Value::INT)
      throw(runtime_error("Execution error - unexpected value on top of stack"));

    return(results.itg);
  }
};


#endif


