#if !defined(LOOS_UNIFORM_WEIGHTS_HPP)
#define LOOS_UNIFORM_WEIGHTS_HPP
#include "Weights.hpp"
#include <Trajectory.hpp>
#include <loos_defs.hpp>

namespace loos {

    class WeightsFromFile : public Weights{
    public:
        uint current_frame;
    private:
        double _total;
        std::string _filename;
        bool _has_list;

    public:
        const double get();
        const double get(const uint index);
        void set(double newWeight);
        void set(double newWeight, const uint index);
        uint size();

        void normalize();
        void accumulate();
        void accumulate(const uint index);
        const double totalWeight();
        const double trajWeight();
        void add_traj(pTraj&  traj);
        const double operator()();
        const double operator()(const uint index);
        void operator()(double newWeight);
        void operator()(double newWeight, const uint index);
        void operator()(std::vector<double>& newWeights);
        
        std::vector<double> weights();

    private:
        uint read_weights(const std::string &filename);
        uint _num_weights;
        pTraj _traj;
        std::vector<double> _weights;
        std::map<std::string, std::string> _weights_files;
        double _totalTraj;

    public:
        WeightsFromFile(const std::string &filename, pTraj& traj ):
                                        current_frame(0),
                                        _total(0.0),
                                        _filename(filename),
                                        _has_list(false)
                                       {
            add_traj(traj);
        };

        WeightsFromFile(const std::string &filename): current_frame(0),
                                              _total(0.0),
                                             _filename(filename),
                                             _has_list(false) {

        };

        WeightsFromFile() : current_frame(0),
                    _total(0.0),
                    _has_list(false)
                    {

        };

        // define virtual destructor inline to ensure vtable gets made correctly.
        virtual ~WeightsFromFile() { }

        uint read_weights_list(const std::string &filename);

    };
}

#endif