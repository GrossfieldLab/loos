#if !defined(LOOS_UNIFORM_WEIGHTS_HPP)
#define LOOS_UNIFORM_WEIGHTS_HPP
#include "Weights.hpp"
#include <Trajectory.hpp>
#include <loos_defs.hpp>

namespace loos {
class UniformWeight : public Weights {
public:
  uint current_frame;
  // were const
  double get();
  double get(const uint index);
  // were not const
  uint size();

  void normalize();
  void accumulate();
  void accumulate(const uint index);
  // were const
  double totalWeight();
  double trajWeight();
  void add_traj(pTraj &traj);
  // were const
  double operator()();
  double operator()(const uint index);
  // this will be built upon request.
  // will be a trajlength vector of 1.0s,
  // so don't ask if you don't need.
  std::vector<double> weights();

private:
  double _frameWeight;
  double _total;
  double _totalTraj;
  std::string _filename;
  std::vector<double> _weights;
  pTraj _traj;
  bool _has_list;

public:
  UniformWeight()
      : current_frame(0), _frameWeight(1.0), _total(0.0), _has_list(false),
        _filename(""){};

  UniformWeight(pTraj &traj)
      : current_frame(0), _frameWeight(1.0), _total(0.0), _filename(""),
        _has_list(false) {
    add_traj(traj);
  };

};
} // namespace loos
#endif