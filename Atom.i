
%include <std_string.i>
%include <std_vector.i>
%include <boost_shared_ptr.i>

%template(IntVector) std::vector<int>;
%shared_ptr(loos::Atom)


%{
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

#include <loos_defs.hpp>
#include "Atom.hpp"

#include <sstream>

%}




namespace loos {

  class Atom;
  typedef boost::shared_ptr<Atom> pAtom;

  class Atom {
  public:
    enum bits {
      nullbit = 0,
      coordsbit = 1,
      bondsbit = coordsbit << 1,
      massbit = bondsbit << 1,
      chargebit = massbit << 1,
      anumbit = chargebit << 1,
      flagbit = anumbit << 1,
      usr1bit = flagbit << 1,
      usr2bit = usr1bit << 1,
      usr3bit = usr2bit << 1
    };

    Atom();
    Atom(const int i, const std::string s, const GCoord& c);
    ~Atom() { }

    int id(void) const;
    void id(const int);
  
    int resid(void) const;
    void resid(const int);

    int atomic_number(void) const;
    void atomic_number(const int);

    std::string name(void) const;
    void name(const std::string);

    std::string altLoc(void) const;
    void altLoc(const std::string);

    std::string chainId(void) const;
    void chainId(const std::string);

    std::string resname(void) const;
    void resname(const std::string);

    std::string segid(void) const;
    void segid(const std::string);

    std::string iCode(void) const;
    void iCode(const std::string);

    std::string PDBelement(void) const;
    void PDBelement(const std::string);

    GCoord& coords(void);
    void coords(const GCoord&);

    double bfactor(void) const;
    void bfactor(const double);

    double occupancy(void) const;
    void occupancy(const double);

    double charge(void) const;
    void charge(const double);

    double mass(void) const;
    void mass(const double);

    std::string recordName(void) const;
    void recordName(const std::string);

    void clearBonds(void);

    void addBond(const pAtom&);
    void addBond(const int);
    void deleteBond(const int);
    void deleteBond(const pAtom&);
    std::vector<int> getBonds(void) const;
    void setBonds(const std::vector<int>& list);
    bool hasBonds(void) const;
    bool isBoundTo(const int);
    bool isBoundTo(const pAtom&);

    bool checkProperty(const bits bitmask);
    void setProperty(const bits bitmask);
    void clearProperty(const bits bitmask);
  };



  %extend Atom {
    char* __str__() {
      static char buf[1024];
      std::ostringstream oss;
      oss << *$self;
      strncpy(buf, oss.str().c_str(), sizeof(buf));
      return(buf);
    }
    char* __repr__() {
      static char buf[1024];
      std::ostringstream oss;
      oss << *$self;
      strncpy(buf, oss.str().c_str(), sizeof(buf));
      return(buf);
    }

    loos::pAtom __copy__() {
      return(loos::pAtom(new loos::Atom(*$self)));
    }

    loos::pAtom __deepcopy(void* p) {
      return(loos::pAtom(new loos::Atom(*$self)));
    }


  };


};

