#if !defined(LOOS_UNIFORM_WEIGHTS_HPP)
#define LOOS_UNIFORM_WEIGHTS_HPP
#include "Weights.hpp"
#include <Trajectory.hpp>
#include <loos_defs.hpp>

namespace loos {

    class WeightsFromFile : protected Weights{
    public:
        uint current_frame;
    private:
        std::string _filename;
        bool _has_list;

    public:
        void add_traj(pTraj&  traj);
        

    private:
        uint read_weights(const std::string &filename);
        std::map<std::string, std::string> _weights_files;

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