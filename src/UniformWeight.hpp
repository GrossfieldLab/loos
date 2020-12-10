#if !defined(LOOS_UNIFORM_WEIGHTS_HPP)
#define LOOS_UNIFORM_WEIGHTS_HPP

#include "Weights.hpp"
#include <Trajectory.hpp>
#include <loos_defs.hpp>

namespace loos {
class UniformWeight : public Weights {
public:
  uint current_frame;
  const double get();
  const double get(const uint index);

  void normalize();
  void accumulate();
  void accumulate(const uint index);
  const double operator()();
  const double operator()(const uint index);
  // this will be built upon request.
  // will be a trajlength vector of 1.0s,
  // so don't ask if you don't need.
  std::vector<double> weights();
private:
  double _frameWeight;

public:
  UniformWeight()
      : current_frame(0), _frameWeight(1.0) {};

  UniformWeight(pTraj &traj)
      : Weights(traj), current_frame(0), _frameWeight(1.0) {};

};
} // namespace loos
#endif