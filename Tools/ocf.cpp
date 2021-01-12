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
  uint max_offset;
  bool group_centroids;
  bool residue_centroids;
  bool com;
  float bondlength;
};
const string indent = "    ";
inline void ocf_at_offset(uint offset, vector<GCoord> &bond_vectors,
                          vector<greal> &outvector) {
  greal accumulated_bvproj = 0;
  greal accumulated_square = 0;
  greal bvproj = 0;
  for (auto i = 0; i < bond_vectors.size() - offset; i++) {
    bvproj = bond_vectors[i].uvdot(bond_vectors[i + offset]);
    accumulated_bvproj += bvproj;
    accumulated_square += bvproj * bvproj;
  }
  outvector[0] = accumulated_bvproj;
  outvector[1] = accumulated_bvproj / (bond_vectors.size() - offset);
  outvector[2] = outvector[1] * outvector[1] -
                 (accumulated_square / (bond_vectors.size() - offset));
}

inline void ocf_at_offset(uint offset, vector<GCoord> &bond_vectors,
                          vector<greal> &outvector, greal weight) {
  greal accumulated_bvproj = 0;
  greal accumulated_square = 0;
  greal bvproj = 0;
  for (auto i = 0; i < bond_vectors.size() - offset; i++) {
    bvproj = bond_vectors[i].uvdot(bond_vectors[i + offset]);
    accumulated_bvproj += bvproj;
    accumulated_square += bvproj * bvproj;
  }
  accumulated_bvproj *= weight;
  outvector[0] = accumulated_bvproj;
  outvector[1] = accumulated_bvproj / (bond_vectors.size() - offset);
  outvector[2] = outvector[1] * outvector[1] -
                 (accumulated_square * weight / (bond_vectors.size() - offset));
}
inline void ag_bond_vectors(AtomicGroup &chain, vector<GCoord> &bond_vectors) {
  for (auto i = 0; i < chain.size() - 1; i++)
    bond_vectors[i] = chain[i]->coords() - chain[i + 1]->coords();
}

inline void centroid_bond_vectors(vector<AtomicGroup> &chain,
                                  vector<GCoord> &bond_vectors) {
  for (auto i = 0; i < chain.size() - 1; i++)
    bond_vectors[i] = chain[i].centroid() - chain[i + 1].centroid();
}

inline void com_bond_vectors(vector<AtomicGroup> &chain,
                             vector<GCoord> &bond_vectors) {
  for (auto i = 0; i < chain.size() - 1; i++)
    bond_vectors[i] = chain[i].centerOfMass() - chain[i + 1].centerOfMass();
}

inline void compute_ocf(uint max_offset, vector<GCoord> &bond_vectors,
                        vector<greal> &ocf_moments_thisframe,
                        vector<greal> &mean_ocfs, vector<greal> &var_ocfs,
                        greal accum_ocfs, greal weight,
                        greal bl_accum) {
  for (auto offset_idx = 0; offset_idx < max_offset; offset_idx++) {
    ocf_at_offset(offset_idx + 1, bond_vectors, ocf_moments_thisframe, weight);
    accum_ocfs += ocf_moments_thisframe[0];
    mean_ocfs[offset_idx] += ocf_moments_thisframe[1];
    var_ocfs[offset_idx] += ocf_moments_thisframe[2];
  }
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
  AtomicGroup scope = selectAtoms(model, sopts->selection);
  pTraj traj = mtopts->trajectory;
  // move unique ptr to Weights into main function ownership for ease of use.
  auto weights = move(wopts->pWeights);
  // Attach trajectory to weights
  weights->addTraj(traj);
  vector<greal> mean_ocfs(topts->max_offset, 0);
  greal bondlength = 0;
  vector<greal> measurement_at_offset(3, 0);
  vector<AtomicGroup> chain;
  if (topts->com) {
    if (topts->group_centroids) {
      chain = scope.splitByMolecule(topts->bond_atom_selection);
      vector<GCoord> bond_vectors(chain.size() - 1, 0);
      for (auto frame_index : mtopts->frameList()) {
        traj->readFrame(frame_index);
        traj->updateGroupCoords(scope);
        // get frame weights; defaults to zero
        const double weight = weights->get();
        weights->accumulate();
        com_bond_vectors(chain, bond_vectors);
        for (auto offset_idx = 0; offset_idx < topts->max_offset;
             offset_idx++) {
          ocf_at_offset(offset_idx + 1, bond_vectors, measurement_at_offset,
                        weight);
          mean_ocfs[offset_idx] += measurement_at_offset[1];
        }
        // compute average bond-length for this frame
        for (auto bond : bond_vectors)
          bondlength += bond.length() / bond_vectors.size();
      }
    } else if (topts->residue_centroids) {
      vector<AtomicGroup> chain =
          selectAtoms(scope, topts->bond_atom_selection).splitByResidue();
      vector<GCoord> bond_vectors(chain.size() - 1, 0);
      for (auto frame_index : mtopts->frameList()) {
        traj->readFrame(frame_index);
        traj->updateGroupCoords(scope);
        // get frame weights; defaults to zero
        const double weight = weights->get();
        weights->accumulate();
        com_bond_vectors(chain, bond_vectors);
        for (auto offset_idx = 0; offset_idx < topts->max_offset;
             offset_idx++) {
          ocf_at_offset(offset_idx + 1, bond_vectors, measurement_at_offset,
                        weight);
          mean_ocfs[offset_idx] += measurement_at_offset[1];
        }
        // compute average bond-length for this frame
        for (auto bond : bond_vectors)
          bondlength += bond.length() * weight / bond_vectors.size();
      }
    }
  } else {
    if (topts->group_centroids) {
      vector<AtomicGroup> chain =
          scope.splitByMolecule(topts->bond_atom_selection);
      vector<GCoord> bond_vectors(chain.size() - 1, 0);
      for (auto frame_index : mtopts->frameList()) {
        traj->readFrame(frame_index);
        traj->updateGroupCoords(scope);
        // get frame weights; defaults to zero
        const double weight = weights->get();
        weights->accumulate();
        centroid_bond_vectors(chain, bond_vectors);
        for (auto offset_idx = 0; offset_idx < topts->max_offset; offset_idx++)
          mean_ocfs[offset_idx] +=
              ocf_at_offset(offset_idx + 1, bond_vectors) * weight;
        // compute average bond-length for this frame
        for (auto bond : bond_vectors)
          bondlength += bond.length() * weight / bond_vectors.size();
      }
    } else if (topts->residue_centroids) {
      vector<AtomicGroup> chain =
          selectAtoms(scope, topts->bond_atom_selection).splitByResidue();
      vector<GCoord> bond_vectors(chain.size() - 1, 0);
      for (auto frame_index : mtopts->frameList()) {
        traj->readFrame(frame_index);
        traj->updateGroupCoords(scope);
        // get frame weights; defaults to zero
        const double weight = weights->get();
        weights->accumulate();
        centroid_bond_vectors(chain, bond_vectors);
        for (auto offset_idx = 0; offset_idx < topts->max_offset; offset_idx++)
          mean_ocfs[offset_idx] +=
              ocf_at_offset(offset_idx + 1, bond_vectors) * weight;
        // compute average bond-length for this frame
        for (auto bond : bond_vectors)
          bondlength += bond.length() * weight / bond_vectors.size();
      }
    } else {
      AtomicGroup chain = selectAtoms(scope, topts->bond_atom_selection);
      vector<GCoord> bond_vectors(chain.size() - 1, 0);
      for (auto frame_index : mtopts->frameList()) {
        traj->readFrame(frame_index);
        traj->updateGroupCoords(scope);
        // get frame weights; defaults to zero
        const double weight = weights->get();
        weights->accumulate();
        ag_bond_vectors(chain, bond_vectors);
        for (auto offset_idx = 0; offset_idx < topts->max_offset; offset_idx++)
          mean_ocfs[offset_idx] +=
              ocf_at_offset(offset_idx + 1, bond_vectors) * weight;
        // compute average bond-length for this frame
        for (auto bond : bond_vectors)
          bondlength += bond.length() * weight / bond_vectors.size();
      }
    }
  }
  cout << "{\n" + indent + "\"mean ocfs\": [\n";
  for (auto i : mean_ocfs)
    cout << indent + indent << i / weights->totalWeight() << ",\n";
  cout << indent + "]\n";
  cout << indent + "\"mean bondlength\": "
       << bondlength / weights->totalWeight() << "\n";
  cout << "}\n";
}
