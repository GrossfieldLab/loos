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
        void addTraj(pTraj &traj);
        

    private:
        uint readWeights(const std::string &filename);
        std::map<std::string, std::string> _weights_files;

    public:
        WeightsFromFile(const std::string &filename, pTraj& traj ):
                                        _filename(filename),
                                        _has_list(false)
                                       {
            addTraj(traj);
        };

        WeightsFromFile(const std::string &filename): 
                                             _filename(filename),
                                             _has_list(false) { };

        WeightsFromFile() : _has_list(false) { };

        // define virtual destructor inline to ensure vtable gets made correctly.
        virtual ~WeightsFromFile() { }

        uint readWeightsList(const std::string &filename);

    };
} // namespace loos

#endif