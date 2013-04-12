
%header %{
#include <loos_defs.hpp>
#include <Coord.hpp>
#include <AtomicGroup.hpp>
#include <PeriodicBox.hpp>
#include <HBondDetector.hpp>
%}

namespace loos {
    class HBondDetector {
    public:
    HBondDetector(const double distance, const double angle, 
                  const AtomicGroup &group);
        
    HBondDetector(const AtomicGroup &group);
    
    HBondDetector();
        
    bool hBonded(const pAtom donor, const pAtom hydrogen, 
                 const pAtom acceptor);
    };
};
