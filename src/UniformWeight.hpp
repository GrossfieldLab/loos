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
  uint size();

  void normalize();
  void accumulate();
  void accumulate(const uint index);
  const double totalWeight();
  const double trajWeight();
  const double operator()();
  const double operator()(const uint index);
  // this will be built upon request.
  // will be a trajlength vector of 1.0s,
  // so don't ask if you don't need.
  std::vector<double> weights();

private:
  double _frameWeight;
  double _total;
  double _totalTraj;
  bool _has_list;
  std::vector<double> _weights;
  pTraj _traj;

public:
  UniformWeight()
      : current_frame(0), _frameWeight(1.0), _total(0.0), _has_list(false) {};

  UniformWeight(pTraj &traj)
      : current_frame(0), _frameWeight(1.0), _total(0.0), _has_list(false), _traj(traj) {};

};
} // namespace loos
#endif