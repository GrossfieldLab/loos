#include <iostream>
#include <loos.hpp>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

const string fullHelpMessage =
    // clang-format off
"XXX";
// clang-format on

class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() {}
  // clang-format off
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("bond-atom-selection,B", po::value<string>(&bond_atom_selection)->
      default_value("name == 'CA' || name == 'P'"),
      "Selection of atoms to compute the OCF across")
      ("max-offset,M", po::value<uint>(&max_offset)->default_value(12),
       "Consider all |i - j| up to this value.")
      ("group-centroids", po::bool_switch(&group_centroids)->default_value(false),
       "If thrown, split bond-atom-selection by molecule and compute BVs between centroids.")
      ("residue-centroids", po::bool_switch(&residue_centroids)->default_value(false),
       "Split bond-atom-selection by residue, then track centroids for bond-vectors.")
      ("infer-connectivity", po::value<float>(&bondlength)->default_value(-1), 
       "Infer connectivity using provided distance for models lacking this. ALERT: "
       "uses hard distance cutoff on first frame of traj to infer connectivity. "
       "Only does this for values greater than zero.")
    ;
  }
  // clang-format on
  string print() const {
    ostringstream oss;
    oss << boost::format("bond_atom_selection=%s,max_offset=%d,group_centroids="
                         "%b,bondlength=%d") %
               bond_atom_selection % max_offset % group_centroids % bondlength;
    return (oss.str());
  }

  bool postcConditions(po::variables_map &map) {
    if (group_centroids && residue_centroids) {
      cerr << "ERROR: --group-centroids and --residue-centroids flags are "
              "mutually exclusive.\n";
      return (false);
    } else
      return (true);
  }
  string bond_atom_selection;
  uint max_offset;
  bool group_centroids;
  bool residue_centroids;
  float bondlength;
};

class OCFAtomic {
private:
  vector<GCoord> bond_vectors;
  vector<AtomicGroup> chain;
  uint max_offset;

public:
  OCFAtomic(vector<AtomicGroup> &input_chain, uint max_offset)
      : bond_vectors(input_chain.size() - 1), max_offset{max_offset} {
    chain = input_chain;
  };
  OCFAtomic(AtomicGroup &input_chain, uint max_offset)
      : bond_vectors(input_chain.size() - 1), max_offset{max_offset} {
    chain = input_chain.splitByMolecule();
  };
  ~OCFAtomic();
  virtual void update_bvs(void) {
    for (auto i = 0; i < chain.size() - 1; i++) {
      bond_vectors[i] = chain[i][0]->coords() - chain[i + 1][0]->coords();
    }
  }
  loos::greal ocf_at_offset(uint offset) {
    loos::greal accumulated_ocf = 0;
    for (auto i = 0; i < bond_vectors.size() - offset; i++) {
      accumulated_ocf += bond_vectors[i].uvdot(bond_vectors[i + offset]);
    }
    return (accumulated_ocf / (bond_vectors.size() - offset));
  }
  vector<loos::greal> compute_all_offsets(void) {
    vector<loos::greal> accumulated_ocf(max_offset);
    update_bvs();
    for (auto offset_idx = 0; offset_idx < max_offset; offset_idx++)
      accumulated_ocf[offset_idx] = ocf_at_offset(offset_idx + 1);
    return accumulated_ocf;
  };
  vector<GCoord> get_bond_vectors() { return (bond_vectors); };
};

OCFAtomic::~OCFAtomic() {}

class OCFCentroids : public OCFAtomic {
private:
  vector<AtomicGroup> chain_groups;
  uint max_offset;
  vector<GCoord> bond_vectors;

public:
  OCFCentroids(vector<AtomicGroup> &input_chains, uint max_offset)
      : bond_vectors(input_chains.size() - 1), max_offset{max_offset} {
    chain_groups = input_chains;
  }
  void update_bvs(void) {
    for (auto i = 0; i < chain_groups.size() - 1; i++)
      bond_vectors[i] =
          chain_groups[i].centroid() - chain_groups[i + 1].centroid();
  }
  ~OCFCentroids();
};

OCFCentroids::~OCFCentroids() {}

class OCFCOM: public OCFAtomic {
private:
  vector<AtomicGroup> chain_groups;
  uint max_offset;
  vector<GCoord> bond_vectors;

public:
  OCFCOM(vector<AtomicGroup> &input_chains, uint max_offset)
      : bond_vectors(input_chains.size() - 1), max_offset{max_offset} {
    chain_groups = input_chains;
  }
  void update_bvs(void) {
    for (auto i = 0; i < chain_groups.size() - 1; i++)
      bond_vectors[i] =
          chain_groups[i].centerOfMass() - chain_groups[i + 1].centerOfMass();
  }
  ~OCFCOM();
};

OCFCOM::~OCFCOM() {};


int main(int argc, char *argv[]) {

  // parse the command line options
  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions *bopts = new opts::BasicOptions(fullHelpMessage);
  opts::BasicSelection *sopts = new opts::BasicSelection("all");
  opts::MultiTrajOptions *mtopts = new opts::MultiTrajOptions;
  opts::WeightsOptions *wopts = new opts::WeightsOptions;
  ToolOptions *topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mtopts).add(wopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  cout << "# " << hdr << "\n";
  // establish system, and subsystems
  AtomicGroup model = mtopts->model;
  if (model.hasBonds()) {
  } else if (topts->bondlength > 0)
    model.findBonds(topts->bondlength);
  else
    throw(LOOSError(
        "Model does not appear to have chemical connectivity, and "
        "infer-connectivity has not been set to a positive value.\n"));
  AtomicGroup scope = selectAtoms(model, sopts->selection);
  pTraj traj = mtopts->trajectory;
  // Attach trajectory to weights
  auto weights = *(wopts->weights);
  vector<AtomicGroup> chains;
  vector<greal> timeseries;
  OCFAtomic *ocf;
  if (topts->group_centroids) {
    chains = scope.splitByMolecule(topts->bond_atom_selection);
     calculator(chains, topts->max_offset);
  }
  if (topts->residue_centroids) {
    chains = scope.splitByResidue(topts->bond_atom_selection);
    ocf_group calculator(chains, topts->max_offset);
  } else {
    // this will make a singleton vector with just the desired AG.

    OCFAtomic calculator(selectAtoms(scope, topts->bond_atom_selection),
                          topts->max_offset);
  }
  greal accum = 0;
  for (auto frame_index : mtopts->frameList()) {
    traj->readFrame(frame_index);
    traj->updateGroupCoords(scope);
    // get frame weights; defaults to zero
    const double weight = weights.get();
    weights.accumulate();
    accum += calculator.compute() * weight;
  }
  cout << accum / weights.totalWeight();
}
