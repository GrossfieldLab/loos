// (c) 2014 Tod D. Romo, Grossfield Lab, URMC

#include <loos.hpp>
#include <algorithm>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


const double prune_factor = 18.0;


// @cond TOOLS_INTERNAL

typedef vector< vector<uint> >     ContactMatrix;


class ToolOptions : public opts::OptionsPackage 
{
public:
  ToolOptions() : threshold(4.0),
		  gnuplot(false)
  {
  }
  

  void addGeneric(po::options_description& o) 
  {
    o.add_options()
      ("threshold,T", po::value<double>(&threshold)->default_value(threshold), "Distance threshold for contact")
      ("gnuplot", po::value<bool>(&gnuplot)->default_value(gnuplot), "Format output for gnuplot");
    
  }

  string print() const 
  {
    ostringstream oss;
    oss << boost::format("threshold=%f,gnuplot=%d") % threshold % gnuplot;
    return(oss.str());
  }
  
  double threshold;
  bool gnuplot;
};


typedef pair<string, string>    SPair;
typedef vector<SPair>         vSPair;

// @endcond

vSPair resmap;

void makeMap() 
{
  resmap.push_back(SPair("GLY", "G"));
  resmap.push_back(SPair("ALA", "A"));
  resmap.push_back(SPair("VAL", "V"));
  resmap.push_back(SPair("LEU", "L"));
  resmap.push_back(SPair("ILE", "I"));
  resmap.push_back(SPair("MET", "M"));
  resmap.push_back(SPair("PHE", "F"));
  resmap.push_back(SPair("TRP", "W"));
  resmap.push_back(SPair("PRO", "P"));
  resmap.push_back(SPair("SER", "S"));
  resmap.push_back(SPair("THR", "T"));
  resmap.push_back(SPair("CYS", "C"));
  resmap.push_back(SPair("TYR", "Y"));
  resmap.push_back(SPair("ASN", "N"));
  resmap.push_back(SPair("GLN", "Q"));
  resmap.push_back(SPair("ASP", "D"));
  resmap.push_back(SPair("GLU", "E"));
  resmap.push_back(SPair("LYS", "K"));
  resmap.push_back(SPair("ARG", "R"));
  resmap.push_back(SPair("HIS", "H"));  
  resmap.push_back(SPair("HSP", "H"));  
}


string lookupAminoAcid(const string& name) 
{
  for (vSPair::const_iterator i = resmap.begin(); i != resmap.end(); ++i)
    if ((*i).first == name)
      return((*i).second);

  return(name);
}


void splotMatrix(const string& hdr, const ContactMatrix& C, const vector<AtomicGroup>& residues, const uint nframes) 
{


  cout << "# " << hdr << endl;
  
  uint n = residues.size();
  vector<string> tags(n);
  for (uint i = 0; i < n; ++i) {
    ostringstream oss;
    oss << lookupAminoAcid(residues[i][0]->resname()) << residues[i][0]->resid();
    tags[i] = (oss.str());
  }
  
  for (uint j=0; j<n; ++j) {
    for (uint i=0; i<n; ++i) {
      cout << j << '\t' << i << '\t' << tags[j] << '\t' << tags[i] << '\t';
      
      double val = (i < j) ? C[j][i] : C[i][j];
      val /= nframes;
      cout << val << endl;
    }
    cout << endl;
  }
}

    
void writeMatrix(const string& hdr, const ContactMatrix& C, const uint nframes) 
{
  uint n = C.size();
  DoubleMatrix M(n, n);
  
  for (uint j=0; j<n; j++) {
    for (uint i=0; i<j; ++i) {
      double d = static_cast<double>(C[j][i]) / nframes;
      M(j, i) = d;
      M(i, j) = d;
    }
    M(j, j) = 1.0;
  }
  

  writeAsciiMatrix(cout, M, hdr);
}



bool findContacts(const AtomicGroup& a, const AtomicGroup& b, const double cutoff) 
{

  if (a.size() > 1 || b.size() > 1) {
    GCoord ac = a.centroid();
    GCoord bc = b.centroid();
    if (ac.distance(bc) > cutoff + prune_factor)
      return(false);

    double cut2 = cutoff * cutoff;

    for (AtomicGroup::const_iterator ai = a.begin(); ai != a.end(); ++ai)
      for (AtomicGroup::const_iterator bi = b.begin(); bi != b.end(); ++bi)
	if ((*ai)->coords().distance2((*bi)->coords()) <= cut2)
	  return(true);
    
    return(false);

  }
  
  return(a[0]->coords().distance(b[0]->coords()) <= cutoff);
}




int main(int argc, char *argv[]) 
{
  
  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts);
  if (! options.parse(argc, argv))
    exit(-1);
  
  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;
  vector<uint> indices = tropts->frameList();

  AtomicGroup subset = selectAtoms(model, sopts->selection);
  double threshold = topts->threshold;
  
  vector<AtomicGroup> residues = subset.splitByResidue();
  uint n = residues.size();
  ContactMatrix C(n, vector<uint>(n, 0));

  cerr << "Processing- ";
  
  
  for (vector<uint>::iterator t = indices.begin(); t != indices.end(); ++t) {
    if ((t - indices.begin()) % 200 == 0)
      cerr << '.';
    
    traj->readFrame(*t);
    traj->updateGroupCoords(model);
    
    for (uint j = 1; j<n; ++j)
      for (uint i=0; i<j; ++i)
	C[j][i] += findContacts(residues[j], residues[i], threshold);
  }

  cerr << "\nDone!\n";
  
  
  if (topts->gnuplot) {
    makeMap();
    splotMatrix(hdr, C, residues, indices.size());
  } else
    writeMatrix(hdr, C, indices.size());
}

