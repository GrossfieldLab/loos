#if !defined(LOOS_WEIGHTS_FROM_FILE_HPP)
#define LOOS_WEIGHTS_FROM_FILE_HPP
#include "Weights.hpp"
#include <Trajectory.hpp>
#include <loos_defs.hpp>

namespace loos {

    class WeightsFromFile : public Weights{
    private:
        std::string _filename;
        bool _has_list;

    public:
        // rely on this to initialize base class.
        void add_traj(pTraj &traj);
        

    private:
        uint read_weights(const std::string &filename);
        std::map<std::string, std::string> _weights_files;

    public:
        WeightsFromFile(const std::string &filename, pTraj& traj ):
                                        _filename(filename),
                                        _has_list(false)
                                       {
            add_traj(traj);
        };

        WeightsFromFile(const std::string &filename): 
                                             _filename(filename),
                                             _has_list(false) { };

        WeightsFromFile() : _has_list(false) { };

        // define virtual destructor inline to ensure vtable gets made correctly.
        virtual ~WeightsFromFile() { }

        uint read_weights_list(const std::string &filename);

    };
} // namespace loos

#endif