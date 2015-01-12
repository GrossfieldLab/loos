// (c) 2014 Tod D. Romo, Grossfield Lab, URMC

#include <loos.hpp>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

// @cond TOOL_INTERNAL
string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "DESCRIPTION\n"
    "\taverager writes out a PDB for the average structure from a trajectory.  If a selection\n"
    "is given (--selection), then the trajectory is first iteratively aligned to an optimal\n"
    "average structure (see aligner).  The '--average' option takes an optional selection that\n"
    "defines what atoms are averaged and written out, otherwise all non-hydrogen and non-solvent\n"
    "atoms are used.  Note that solvent is selected by a segid of either 'BULK' or 'SOLVENT'.\n"
    "If your system uses a different identifier, you will want to explicitly give a selection\n"
    "for the --average option\n"
    "EXAMPLES\n"
    "\n"
    "\n"
    "SEE ALSO\n";

  return(msg);
}


class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("id,I", po::value<bool>(&id)->default_value(false), "Atom ID")
      ("name,T", po::value<bool>(&name)->default_value(false), "Atom name")
      ("resid,R", po::value<bool>(&resid)->default_value(false), "Residue ID")
      ("resname,N", po::value<bool>(&resname)->default_value(false), "Residue name")
      ("segid,S", po::value<bool>(&segid)->default_value(true), "Segid or Segname")
      ("all,A", "Use all metadata");
  }


  bool postConditions(po::variables_map& map) {
    if (map.count("all")) {
      id = true;
      name = true;
      resname = true;
      resid = true;
      segid = true;
    }
    return(true);
  }

  string print() const {
    ostringstream oss;

    oss << boost::format("id=%d,name=%d,resid=%d,resname=%d,segid=%d") %
      id % name % resname % resid % segid;
    return(oss.str());
  }

  bool id, name, resname, resid, segid;
};

class Tracker {
public:
  Tracker() : _count(0) { }
  

  virtual ~Tracker() {}

  virtual void add(const pAtom& a) =0;
  virtual string report() const =0;

protected:
  string formatLabel(const string& label, const uint width=40) const {
    int remainder = width - 4 - label.size();
    if (remainder < 0)
      return(label);

    string out = "== " + label + " ";
    while (remainder-- > 0)
      out += '=';

    return(out);
  }

  uint _count;
};


class IntTracker : public Tracker {
public:
  IntTracker() : _min(numeric_limits<int>::max()),_max(numeric_limits<int>::min()) {}
  IntTracker(const string& label) : _min(numeric_limits<int>::max()),_max(numeric_limits<int>::min()),
				    _label(label) {}

  virtual void add(const pAtom& a) {
    int i = retrieveAttribute(a);
    if (i < _min)
      _min = i;
    if (i > _max)
      _max = i;
  }


  virtual string report() const {
    ostringstream oss;
    oss << boost::format("%s\n  %10s = %d\n  %10s = %d\n  %10s = %d\n")
      % formatLabel(_label) 
      % "min" % _min
      % "max" % _max
      % "count" % _count;
    return(oss.str());
  }

private:
  virtual int retrieveAttribute(const pAtom& a) =0;

private:
  int _min, _max;
  string _label;
};



class AtomidTracker : public IntTracker {
public:
  AtomidTracker() : IntTracker("Atom ID") { }
private:
  virtual int retrieveAttribute(const pAtom& a) { 
    ++_count; 
    return(a->id());
  }
};


class ResidTracker : public IntTracker {
public:
  ResidTracker() : IntTracker("Residue ID"), _last_resid(-999999) { }
  
private:
  virtual int retrieveAttribute(const pAtom& a) { 
    int i = a->resid();
    if (i != _last_resid) {
      ++_count;
      _last_resid = i;
    }

    return(a->resid());
  }

  int _last_resid;
};


// --------------------------------------------------------------

class UniqueStringTracker : public Tracker {
public:
  UniqueStringTracker(const string& label) : _label(label) {}

  virtual void add(const pAtom& a) {
    ++_count;
    string what = retrieveAttribute(a);
    map<string, uint>::iterator i = _vals.find(what);
    if (i != _vals.end())
      _vals[what] += 1;
    else
      _vals[what] = 1;
  }

  virtual string report() const {
    ostringstream oss;
    oss << formatLabel(_label) << "\n";
    oss << boost::format("* Number of atoms: %d\n") % _count;
    oss << boost::format("* Number of unique values: %d\n") % _vals.size();
    for (map<string, uint>::const_iterator i = _vals.begin(); i != _vals.end(); ++i)
      oss << boost::format("  %10s = %-8d (%.2f %%)\n") % i->first % i->second
	% (100.0 * i->second / _count);
    return(oss.str());
  }
 
private:
  virtual string retrieveAttribute(const pAtom& a) = 0;


private:
  string _label;
  map<string, uint> _vals;

};

class ResnameTracker : public UniqueStringTracker {
public:
  ResnameTracker() : UniqueStringTracker("Residue Name") {}
private:
  virtual string retrieveAttribute(const pAtom& a) { return(a->resname()); }
};

class NameTracker : public UniqueStringTracker {
public:
  NameTracker() : UniqueStringTracker("Atom Name") {}
private:
  virtual string retrieveAttribute(const pAtom& a) { return(a->name()); }
};


class SegidTracker : public UniqueStringTracker {
public:
  SegidTracker() : UniqueStringTracker("Segment Name") {}
private:
  virtual string retrieveAttribute(const pAtom& a) { return(a->segid()); }
};

// --------------------------------------------------------------

class MetaTracker {
public:

  void addTracker(Tracker* t) { _trackers.push_back(t); }
  void processAtom(const pAtom& a) {
    for (vector<Tracker*>::iterator i = _trackers.begin(); i != _trackers.end(); ++i)
      (*i)->add(a);
  }

  string report() const {
    ostringstream oss;
    for (vector<Tracker*>::const_iterator i = _trackers.begin(); i != _trackers.end(); ++i)
      oss << (*i)->report() << '\n';
     
    return(oss.str());
  }
  
private:
  vector<Tracker*> _trackers;

};


// @endcond

int main(int argc, char* argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::RequiredArguments* ropts = new opts::RequiredArguments;
  ropts->addArgument("model", "Model Filename");
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(topts).add(ropts);
  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup model = createSystem(ropts->value("model"));
  AtomicGroup subset = selectAtoms(model, sopts->selection);

  MetaTracker mt;
  if (topts->id)
    mt.addTracker(new AtomidTracker);
  if (topts->name)
    mt.addTracker(new NameTracker);
  if (topts->resname)
    mt.addTracker(new ResnameTracker);
  if (topts->resid)
    mt.addTracker(new ResidTracker);
  if (topts->segid)
    mt.addTracker(new SegidTracker);

  for (AtomicGroup::iterator i = subset.begin(); i != subset.end(); ++i)
    mt.processAtom(*i);

  cout << mt.report();
}
