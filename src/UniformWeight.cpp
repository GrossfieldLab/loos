#include "UniformWeight.hpp"

namespace loos {

// Add trajectory to the instance of UniformWeight
void UniformWeight::add_traj(pTraj &traj) {
  _traj = traj;
  _totalTraj = 0.0;
}

// Accumulate the 'weight' used so far
void UniformWeight::accumulate() {
  _total += _frameWeight;
  _totalTraj += _frameWeight;
}

void UniformWeight::accumulate(const uint index) {
  _total += _frameWeight;
  _totalTraj += _frameWeight;
}

// Normalize weights using length of traj
void UniformWeight::normalize() {
  _frameWeight /= _traj->nframes();
  if (!_weights.empty()) {
    for (uint i = 0; i < _weights.size(); i++)
      _weights[i] = _frameWeight;
  }
}

// return totalWeight used
double UniformWeight::totalWeight() { return _total; }

// return the weight of current trajectory
double UniformWeight::trajWeight() { return _totalTraj; }

// return a vector with the number of weights in it
std::vector<double> UniformWeight::weights() {
  std::vector<double> uniform(_frameWeight, _traj->nframes());
  _weights = uniform;
  return _weights;
}

// return the weighting of the current frame
double UniformWeight::get() {
  current_frame = _traj->currentFrame();
  return _frameWeight;
}

// return the weight of the index
double UniformWeight::get(const uint index) { return _frameWeight; }

// calling nomenclature wraps get
double UniformWeight::operator()() { return get(); }

double UniformWeight::operator()(const uint index) { return get(index); }

// give number of frames, which is equivalent to _num_weights, as size
uint UniformWeight::size() { return _traj->nframes(); }

} // namespace loos