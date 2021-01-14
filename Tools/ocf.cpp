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
      ("max-offset,M", po::value<int>(&max_offset)->default_value(-1),
       "Consider all |i - j| up to this value. -1 considers all possible offsets.")
      ("group-centroids", po::bool_switch(&group_centroids)->default_value(false),
       "If thrown, split bond-atom-selection by molecule and compute BVs between centroids.")
      ("residue-centroids", po::bool_switch(&residue_centroids)->default_value(false),
       "Split bond-atom-selection by residue, then track centroids for bond-vectors.")
      ("center-of-mass,c", po::bool_switch(&com)->default_value(false),
       "Instead of using centroids, use centers of mass for groups/residues.")
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
                         "%b,bondlength=%d,residue_centroids%b,com=%b") %
               bond_atom_selection % max_offset % group_centroids % bondlength %
               residue_centroids % com;
    return (oss.str());
  }

  bool postConditions(po::variables_map &map) {
    if (group_centroids && residue_centroids) {
      cerr << "ERROR: --group-centroids and --residue-centroids flags are "
              "mutually exclusive.\n";
      return (false);
    } else if (com && !(group_centroids || residue_centroids)) {
      cerr << "ERROR: --center-of-mass must be used with --group-centroids or"
              "--residue-centroids.\n";
      return (false);
    } else
      return (true);
  }
  string bond_atom_selection;
  int max_offset;
  bool group_centroids;
  bool residue_centroids;
  bool com;
  float bondlength;
};
const string indent = "    ";

inline void ag_bond_vectors(AtomicGroup &chain, vector<GCoord> &bond_vectors) {
  for (uint i = 0; i < chain.size() - 1; i++)
    bond_vectors[i] = chain[i]->coords() - chain[i + 1]->coords();
}

inline void centroid_bond_vectors(vector<AtomicGroup> &chain,
                                  vector<GCoord> &bond_vectors) {
  for (uint i = 0; i < chain.size() - 1; i++)
    bond_vectors[i] = chain[i].centroid() - chain[i + 1].centroid();
}

inline void com_bond_vectors(vector<AtomicGroup> &chain,
                             vector<GCoord> &bond_vectors) {
  for (uint i = 0; i < chain.size() - 1; i++)
    bond_vectors[i] = chain[i].centerOfMass() - chain[i + 1].centerOfMass();
}
// this is the work to be done inside the traj loop, that is, per-frame.
inline void compute_ocf_bondlength(uint max_offset,
                                   vector<GCoord> &bond_vectors,
                                   greal &accum_ocfs, vector<greal> &mean_ocfs,
                                   vector<greal> &var_ocfs,
                                   vector<greal> &accum_sq_mean, greal weight,
                                   greal &bl_accum) {
  for (uint offset_idx = 0; offset_idx < max_offset; offset_idx++) {
    uint offset = offset_idx + 1;
    greal accumulated_bvproj = 0;
    greal accumulated_square = 0;
    greal bvproj = 0;
    for (uint i = 0; i < bond_vectors.size() - offset; i++) {
      bvproj = bond_vectors[i].uvdot(bond_vectors[i + offset]);
      accumulated_bvproj += bvproj;
      accumulated_square += bvproj * bvproj;
    }
    accum_ocfs += accumulated_bvproj * weight;
    greal mean_ocf_atoffset =
        accumulated_bvproj / (bond_vectors.size() - offset) * weight;
    mean_ocfs[offset_idx] += mean_ocf_atoffset;
    var_ocfs[offset_idx] +=
        accumulated_square * weight / (bond_vectors.size() - offset) -
        mean_ocf_atoffset * mean_ocf_atoffset;
    accum_sq_mean[offset_idx] += mean_ocf_atoffset * mean_ocf_atoffset; 
  }
  for (auto bond : bond_vectors)
    bl_accum += bond.length() * weight;
}

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
  if (topts->max_offset == 0)
    throw(
        LOOSError("You asked for an offset of zero, which is not possible.\n"));
  AtomicGroup scope = selectAtoms(model, sopts->selection);
  pTraj traj = mtopts->trajectory;
  // move unique ptr to Weights into main function ownership for ease of use.
  auto weights = move(wopts->pWeights);
  // Attach trajectory to weights
  weights->addTraj(traj);
  // initialize max offset at top level, define either with user input
  // or as a function of chain, below.
  uint max_offset;
  vector<greal> mean_ocfs;
  vector<greal> var_ocfs;
  vector<greal> accum_sq_mean;
  greal accum_ocf = 0;
  greal bondlength = 0;
  vector<AtomicGroup> chain;
  if (topts->com) {
    if (topts->group_centroids) {
      chain = scope.splitByMolecule(topts->bond_atom_selection);
      vector<GCoord> bond_vectors(chain.size() - 1, 0);
      if (topts->max_offset > 0)
        max_offset = topts->max_offset;
      else if (topts->max_offset < 0)
        max_offset = bond_vectors.size() - 1;
      mean_ocfs.resize(max_offset, 0);
      var_ocfs.resize(max_offset, 0);
      accum_sq_mean.resize(max_offset, 0);
      for (auto frame_index : mtopts->frameList()) {
        traj->readFrame(frame_index);
        traj->updateGroupCoords(scope);
        // get frame weights; defaults to zero
        const double weight = weights->get();
        weights->accumulate();
        com_bond_vectors(chain, bond_vectors);
        compute_ocf_bondlength(max_offset, bond_vectors, accum_ocf, mean_ocfs,
                               var_ocfs, accum_sq_mean, weight, bondlength);
      }
      bondlength /= bond_vectors.size() * mtopts->frameList().size();
    } else if (topts->residue_centroids) {
      vector<AtomicGroup> chain =
          selectAtoms(scope, topts->bond_atom_selection).splitByResidue();
      if (topts->max_offset > 0)
        max_offset = topts->max_offset;
      else if (topts->max_offset < 0)
        max_offset = 0;
      mean_ocfs.resize(max_offset, 0);
      var_ocfs.resize(max_offset, 0);
      accum_sq_mean.resize(max_offset, 0);
      vector<GCoord> bond_vectors(chain.size() - 1, 0);
      for (auto frame_index : mtopts->frameList()) {
        traj->readFrame(frame_index);
        traj->updateGroupCoords(scope);
        // get frame weights; defaults to zero
        const double weight = weights->get();
        weights->accumulate();
        com_bond_vectors(chain, bond_vectors);
        compute_ocf_bondlength(max_offset, bond_vectors, accum_ocf, mean_ocfs,
                               var_ocfs, accum_sq_mean, weight, bondlength);
      }
      bondlength /= bond_vectors.size() * mtopts->frameList().size();
    }
  } else {
    if (topts->group_centroids) {
      vector<AtomicGroup> chain =
          scope.splitByMolecule(topts->bond_atom_selection);
      if (topts->max_offset > 0)
        max_offset = topts->max_offset;
      else if (topts->max_offset < 0)
        max_offset = chain.size() - 2;
      mean_ocfs.resize(max_offset, 0);
      var_ocfs.resize(max_offset, 0);
      accum_sq_mean.resize(max_offset, 0);
      vector<GCoord> bond_vectors(chain.size() - 1, 0);
      for (auto frame_index : mtopts->frameList()) {
        traj->readFrame(frame_index);
        traj->updateGroupCoords(scope);
        // get frame weights; defaults to zero
        const double weight = weights->get();
        weights->accumulate();
        centroid_bond_vectors(chain, bond_vectors);
        compute_ocf_bondlength(max_offset, bond_vectors, accum_ocf, mean_ocfs,
                               var_ocfs, accum_sq_mean, weight, bondlength);
      }
      bondlength /= bond_vectors.size() * mtopts->frameList().size();
    } else if (topts->residue_centroids) {
      vector<AtomicGroup> chain =
          selectAtoms(scope, topts->bond_atom_selection).splitByResidue();
      if (topts->max_offset > 0)
        max_offset = topts->max_offset;
      else if (topts->max_offset < 0)
        max_offset = chain.size() - 2;
      mean_ocfs.resize(max_offset, 0);
      var_ocfs.resize(max_offset, 0);
      accum_sq_mean.resize(max_offset, 0);
      vector<GCoord> bond_vectors(chain.size() - 1, 0);
      for (auto frame_index : mtopts->frameList()) {
        traj->readFrame(frame_index);
        traj->updateGroupCoords(scope);
        // get frame weights; defaults to zero
        const double weight = weights->get();
        weights->accumulate();
        centroid_bond_vectors(chain, bond_vectors);
        compute_ocf_bondlength(max_offset, bond_vectors, accum_ocf, mean_ocfs,
                               var_ocfs, accum_sq_mean, weight, bondlength);
      }
      bondlength /= bond_vectors.size() * mtopts->frameList().size();
    } else {
      AtomicGroup chain = selectAtoms(scope, topts->bond_atom_selection);
      vector<GCoord> bond_vectors(chain.size() - 1, 0);
      if (topts->max_offset > 0)
        max_offset = topts->max_offset;
      else if (topts->max_offset < 0)
        max_offset = chain.size() - 2;
      mean_ocfs.resize(max_offset, 0);
      var_ocfs.resize(max_offset, 0);
      accum_sq_mean.resize(max_offset, 0);
      for (auto frame_index : mtopts->frameList()) {
        traj->readFrame(frame_index);
        traj->updateGroupCoords(scope);
        // get frame weights; defaults to zero
        const double weight = weights->get();
        weights->accumulate();
        ag_bond_vectors(chain, bond_vectors);
        compute_ocf_bondlength(max_offset, bond_vectors, accum_ocf, mean_ocfs,
                               var_ocfs, accum_sq_mean, weight, bondlength);
      }
      bondlength /= bond_vectors.size();
    }
  }
  cout << "{\n" + indent + "\"mean ocfs\": [\n";
  for (auto i : mean_ocfs)
    cout << indent + indent << i / weights->totalWeight() << ",\n";
  cout << indent + "],\n";
  cout << indent + "\"variance of means\": [\n";
  greal mean_ocf;
  for (uint i = 0; i < mean_ocfs.size(); i++) {
    mean_ocf = mean_ocfs[i] / weights->totalWeight();
    cout << indent + indent
         << accum_sq_mean[i] / weights->totalWeight() - mean_ocf * mean_ocf
         << ",\n";
  }
  cout << indent + "],\n";
  cout << indent + "\"mean variances\": [\n";
  for (auto i : var_ocfs)
    cout << indent + indent << i / weights->totalWeight() << ",\n";
  cout << indent + "],\n";
  cout << indent + "\"mean projections summed\": "
       << accum_ocf / weights->totalWeight() << ",\n";
  cout << indent + "\"mean bondlength\": "
       << bondlength / weights->totalWeight() << "\n";
  cout << "}\n";
}
