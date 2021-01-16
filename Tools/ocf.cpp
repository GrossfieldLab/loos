#include <iostream>
#include <loos.hpp>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

const string fullHelpMessage =
    // clang-format off
"SYNOPSIS \n"
" \n"
"This tool is designed to compute the Orientational Correlation Function, in the\n"
" style of its use in polymer chemistry contexts such as in Ullner, M. & \n"
"Woodward, C. E. Orientational Correlation Function and Persistence Lengths of \n"
"Flexible Polyelectrolytes. Macromolecules 35, 1437–1445 (2002) and more \n"
"specifically as it was used in Plumridge, A., Andresen, K. & Pollack, L. \n"
"Visualizing Disordered Single-Stranded RNA: Connecting Sequence, Structure, and\n"
" Electrostatics. J. Am. Chem. Soc. 142, 109–119 (2020). The user may specify \n"
"abstracted 'bond vectors' between links in the polymer chain using single atom \n"
"selectors or group selectors.  \n"
" \n"
"DESCRIPTION \n"
" \n"
"This tool uses the definition of the orientational correlation function from \n"
"Ullner & Woodward to estimate how correlated links in a polymer chain are. This\n"
" is done by looking at the normalized projection of the i bond-vector in the \n"
"chain onto the i+o bond-vector in the chain, for all offsets o between 1 and a \n"
"max offset specified by the user (default is -1, or all possible). Each bond \n"
"vector is defined as being a link between a point on a certain residue and a \n"
"point on a neighboring residue. These could literally be a bond vector, if the \n"
"points are atoms bonded to one another, or it could be a 'coarse grained' \n"
"linkage between two monomers in the chain. For example, in the CA default for \n"
"proteins, each residue is being treated as a link, with the link position at \n"
"the alpha carbon; in such a coarsening of the polypeptide chain the vector \n"
"between CAs becomes the chain bond vector.  \n"
" \n"
"Thus, in the CA example, it would be the projection of the vector between CA_i \n"
"and CA_i+1, v_i, onto the vector between CA_{i+o} and CA_{i+o+1}. These \n"
"projections are averaged across all possible i for each o requested, then \n"
"reported as a list of correlations as a function of offset. It is also possible\n"
" for anticorrelations to be exhibited by this function--for example, a pretty \n"
"solid antiparallel beta sheet would likely produce vectors that are pointed in \n"
"opposite directions but are nearly coplanar. \n"
" \n"
"The notion here is that stiffer chains have a persistence of orientation, which\n"
" is quantified by the projection of these vectors. Thus, a 'length' is also \n"
"defined; it is the average length of the bond-vectors, multiplied by the \n"
"average correlation between bond vectors summed over all pairs.  \n"
" \n"
"The tool writes the results of the requested analysis to stdout as JSON. The \n"
"JSON has the following tags: \"mean ocfs\", \"variance of means\", \"mean \n"
"variances\", \"mean projections summed\", and \"mean bondlength\". The first three \n"
"are all arrays with lengths equal to the number of offsets analyzed. The \"mean \n"
"ocfs\" are the normalized projection vectors averaged over all pairs of bond \n"
"vectors with a given offset, and then across each frame analyzed. The \"variance\n"
" of means\" is the variance in each such mean across the whole trajectory. The \n"
"\"mean variances\" are the variances at each offset, averaged across all analyzed\n"
" frames. Finally, the \"mean projections summed\" and \"mean bondlength\" when \n"
"multiplied together should correspond to l_OCF as given in Plumridge et al. \n"
" \n"
"The idea behind reporting the mean variances and the variances of the mean is \n"
"that the one reports on the variability of each chain projection across the \n"
"trajectory, whereas the other reports on how variable the projections that are \n"
"associated with each offset are within a given frame, on average. If the \n"
"variability within an offset is high on average, but the variability of that \n"
"offset is low across the trajectory, it implies a strong conformational \n"
"preference, or static disorder/glassy behavior. This may mean the chain is not \n"
"really sampling different conformations, even if it is exhibiting correlation \n"
"die-off as a function of offset. \n"
" \n"
"EXAMPLES \n"
" \n"
"ocf model traj > ocf_traj.json \n"
" \n"
"This will look for either alpha carbons or phosphorus atoms within the entire \n"
"model, then use their ordering in the model to draw vectors between each such \n"
"atom and the next one in the chain. It will compute the ocf on each frame in \n"
"the trajectory. It will do this for all possible offsets. \n"
" \n"
"ocf --bond-atom-selection 'name ~= \"C*\'\"' --center-of-mass --residue-centroids\n"
" \ \n"
"model traj > ocf_sugar_carbons_com.json \n"
" \n"
"This will use the centers of mass of atoms matching the regex 'C*\'' (any \n"
"carbon with a single quote at the end, which hopefully amounts to sugar \n"
"carbons) as the points between which to draw bond-vectors for each link in the \n"
"chain. It will then proceed to compute the OCF as normal for these. To use the \n"
"centroids instead of the COM, elide the --center-of-mass flag.  \n"
" \n"
"ocf --selection 'resid < 31' --bond-atom-selection 'name ~= \"C*\'\"' --center-\n"
"of-mass --residue-centroids \ \n"
"model traj > ocf_sugar_carbons_com.json \n"
" \n"
"Like the above, but only look for the bond atom selection within the first 30 \n"
"residues in the model. \n"
" \n"
"POTENTIAL COMPLICATIONS \n"
" \n"
"Be careful with selection strings; results that are only subtly wrong could \n"
"emerge from a string that grabs atoms or groups you're not expecting. While \n"
"this is always a good thing to be careful about when analyzing trajectories, \n"
"the peril here (because the selections are being split internally across either\n"
" residues or contiguous sections within the bond atom selection) seems great \n"
"indeed. \n"
" \n"
"The '--group-centroids' flag shouldn't be used unless you're after treating a \n"
"collection of atoms that is trans-residue according to how your model defines \n"
"residues. If you do need this functionality, make sure your model has \n"
"connectivity, or find a way to add it. Using the '--infer-connectivity' flag to\n"
" do this is applying a simple distance cutoff to decide where the chemical \n"
"bonds are in your system from the first frame, which should be treated with \n"
"caution.\n"
;
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
  vector<AtomicGroup> chain;
  // Define pointer to the function that will update the bond sites -- will be
  // either com_bond_vectors or centroid_bond_vectors, depending on user input
  void (*bv_getter)(vector<AtomicGroup> &, vector<GCoord> &);
  // determine what points to use for the centers of each link in the chain
  // based on user input
  if (topts->com) {
    bv_getter = &com_bond_vectors;
    if (topts->group_centroids) {
      chain = scope.splitByMolecule(topts->bond_atom_selection);
    } else if (topts->residue_centroids) {
      chain = selectAtoms(scope, topts->bond_atom_selection).splitByResidue();
    }
  } else {
    bv_getter = &centroid_bond_vectors;
    if (topts->group_centroids) {
      chain = scope.splitByMolecule(topts->bond_atom_selection);
    } else if (topts->residue_centroids) {
      chain = selectAtoms(scope, topts->bond_atom_selection).splitByResidue();
    } else {
      chain = selectAtoms(scope, topts->bond_atom_selection).splitByMolecule();
    }
  }
  // Now figure out how many bond vectors and offsets we're tracking
  vector<GCoord> bond_vectors(chain.size() - 1, 0);
  vector<greal> mean_ocfs(chain.size() - 2, 0);
  vector<greal> var_ocfs(chain.size() - 2, 0);
  vector<greal> accum_sq_mean(chain.size() - 2, 0);
  greal accum_ocf = 0;
  greal bondlength = 0;
  uint max_offset;
  if (topts->max_offset > 0)
    max_offset = topts->max_offset;
  else if (topts->max_offset < 0)
    max_offset = bond_vectors.size() - 1;
  // loop over trajectory
  for (auto frame_index : mtopts->frameList()) {
    traj->readFrame(frame_index);
    traj->updateGroupCoords(scope);
    // get frame weights; defaults to zero
    const double weight = weights->get();
    weights->accumulate();
    bv_getter(chain, bond_vectors);
    compute_ocf_bondlength(max_offset, bond_vectors, accum_ocf, mean_ocfs,
                           var_ocfs, accum_sq_mean, weight, bondlength);
  }
  bondlength /= bond_vectors.size();
  // create the JSON report, written to stdout.
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
