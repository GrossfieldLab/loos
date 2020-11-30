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

inline greal by_atom_ocf(vector<AtomicGroup> &ag_singleton, uint max_offset) {
  greal accumulated_ocf = 0;
  for (auto offset = 1; offset == max_offset; offset++) {
    accumulated_ocf += ag_singleton[0].ocf(offset);
  }
  return (accumulated_ocf / (max_offset - 1));
}

inline greal by_group_centroid_ocf(vector<AtomicGroup> &ag_list,
                                   uint max_offset) {
  greal accumulated_ocf = 0;
  vector<GCoord> bvs;
  for (auto i = 0; i < ag_list.size() - 1; i++)
    bvs.push_back(ag_list[i].centroid() - ag_list[i + 1].centroid());
  for (auto offset = 1; offset < max_offset + 1; offset++) {
    for (auto i = 0; i < bvs.size() - offset; i++)
      accumulated_ocf = bvs[i].dot(bvs[i + offset]) /
                        (bvs[i].length() * bvs[i + offset].length());
  }
  return (accumulated_ocf / (bvs.size() * max_offset));
}

inline greal mean_bond_length(vector<AtomicGroup>)