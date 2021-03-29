#ifndef __FORMFACTORSET_HPP
#define __FORMFACTORSET_HPP

/*
 *   Store a collection for FormFactor objects for the nuclei have them for
 */

#include <map>
#include <exceptions.hpp>
#include <FormFactor.hpp>

namespace loos {

    class FormFactorSet
    {

    public:
        FormFactorSet()
            {
            setup();
            }

        double get(uint i, double q) {
            std::map<uint,FormFactor>::iterator it = _map.find(i);
            if (it == _map.end()) {
                throw(LOOSError(" unsupported atomic number in scattering calculation"));
            }
            double val = _map[i].compute(q);
            return val;
        }

    private:
        std::map<uint,FormFactor> _map;
        void setup();
    };
}



#endif
