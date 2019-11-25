#include "UniformWeight.hpp"

namespace loos {
    // return a vector with the number of weights in it
    std::vector<double> UniformWeight::weights() {
        std::vector<double> uniform(1.0, _num_weights);
        _weights = uniform;
        return _weights;
    }

    // return the weighting of the current frame
    const double UniformWeight::get() {
        return 1.0;
    }

    // return the weight of the current index
    const double UniformWeight::get() {
        return 1.0;
    }

    // calling nomenclature wraps get
    const double Weights::operator()(){
        return get();
    }

    // calling nomenclature wraps get
    const double Weights::operator(const uint index)(){
        return get(index);
    }

    // give number of frames, which is equivalent to _num_weights, as size
    uint UniformWeights::size() {
        return _num_weights;
    }

    // Return the total weight used so far (from accumulate) 
}