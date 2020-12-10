#include "UniformWeight.hpp"

namespace loos {

//! Accumulate the 'weight' used so far
void UniformWeight::accumulate() {
  _total += _frameWeight;
  _totalTraj += _frameWeight;
}
//! accumulate the 'weight' used so far.
//! Do nothing with index, because weights are constant.
void UniformWeight::accumulate(const uint index) {
  _total += _frameWeight;
  _totalTraj += _frameWeight;
}

//! Normalize weights using length of traj
void UniformWeight::normalize() {
  _frameWeight /= _traj->nframes();
  if (!_weights.empty()) {
    for (uint i = 0; i < _weights.size(); i++)
      _weights[i] = _frameWeight;
  }
}

//! return totalWeight used
const double UniformWeight::totalWeight() { return _total; }

//! return the weight of current trajectory
const double UniformWeight::trajWeight() { return _totalTraj; }

//! return a vector with the number of weights in it
std::vector<double> UniformWeight::weights() {
  std::vector<double> uniform(_frameWeight, _traj->nframes());
  _weights = uniform;
  return _weights;
}

//! return the weighting of the current frame
const double UniformWeight::get() {
  current_frame = _traj->currentFrame();
  return _frameWeight;
}

//! return the weight of the index
const double UniformWeight::get(const uint index) { return _frameWeight; }

//! calling nomenclature wraps get
const double UniformWeight::operator()() { return get(); }

//! ignore the index, frame weights are constant.
const double UniformWeight::operator()(const uint index) { return get(); }

} // namespace loos