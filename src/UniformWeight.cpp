#include "UniformWeight.hpp"

namespace loos
{

// Accumulate the 'weight' used so far

// Normalize weights using length of traj
void UniformWeight::normalize()
{
    _frameWeight /= _traj->nframes();
    if (!_weights.empty())
    {
        for (uint i = 0; i < _weights.size(); i++)
            _weights[i] = _frameWeight;
    }
}

// return a vector with the number of weights in it
std::vector<double> UniformWeight::weights()
{
    std::vector<double> uniform(_frameWeight, _traj->nframes());
    _weights = uniform;
    return _weights;
}

// return the weighting of the current frame
const double UniformWeight::get()
{
    current_frame = _traj->currentFrame();
    return _frameWeight;
}

// return the weight of the current index
const double UniformWeight::get(const uint index)
{
    return _frameWeight;
}

// calling nomenclature wraps get
const double UniformWeight::operator()()
{
    return get();
}

// calling nomenclature wraps get
const double UniformWeight::operator()(const uint index)
{
    return get(index);
}

// give number of frames, which is equivalent to _num_weights, as size
uint UniformWeight::size()
{
    return _traj->nframes();
}
} // namespace loos