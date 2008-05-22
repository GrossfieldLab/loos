/*
  Selectors.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  A selector library (of sorts...)
*/




#if !defined(SELECTORS_HPP)
#define SELECTORS_HPP


#include <AtomicGroup.hpp>
#include <Kernel.hpp>



struct CAlphaSelector : public AtomSelector {
  bool operator()(const pAtom& pa) const {
    return(pa->name() == "CA");
  }
};


struct SolventSelector : public AtomSelector {
  bool operator()(const pAtom& pa) const {
    return(pa->segid() == "SOLV" || pa->segid() == "BULK");
  }
};

struct MainChainSelector : public AtomSelector {
  bool operator()(const pAtom& pa) const {
    string s = pa->name();
    return(s == "C" || s == "CA" || s == "O" || s == "N");
  }
};


struct SegidSelector : public AtomSelector {
  SegidSelector(const string s) : str(s) { }
  bool operator()(const pAtom& pa) const {
    return(pa->segid() == str);
  }

  string str;
};


struct NotSelector : public AtomSelector {
  NotSelector(const AtomSelector& s) : sel(s) { }
  bool operator()(const pAtom& pa) const {
    return(!(sel(pa)));
  }

  const AtomSelector& sel;
};




class CompiledSelector : public AtomSelector {
public:
  CompiledSelector(loos::Kernel& k) : krnl(k) { }

  bool operator()(const pAtom& pa) const {
    krnl.execute(pa);
    if (krnl.stack().size() != 1) {
      throw(runtime_error("Execution error - unexpected values on stack"));
    }

    loos::Value results = krnl.stack().pop();
    if (results.type != loos::Value::INT)
      throw(runtime_error("Execution error - unexpected value on top of stack"));

    return(results.itg);
  }


private:
  loos::Kernel& krnl;

};




#endif

