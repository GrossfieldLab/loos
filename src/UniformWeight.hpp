#if !defined(LOOS_UNIFORM_WEIGHTS_HPP)
#define LOOS_UNIFORM_WEIGHTS_HPP
#include "Weights.hpp"

namespace loos {
    class UniformWeight : public Weights {
    public:
        const double get();
        const double get(const uint index);
        uint size();

        void normalize();
        void accumulate();
        void accumulate(const uint index);
        const double totalWeight();
        const double trajWeight();
        void add_traj(ptraj& traj);
        std::vector<double> weights();
    };
    private:
        double _total;
        std::string _filename;
        std::vector<double> _weights;
        pTraj _traj;
        bool _has_list;
    
    public:
    UniformWeight() : current_frame(0),
                      _total(0.0),
                      _has_list(false),
                      _filename("") {};
    ~UniformWeight() { };
}