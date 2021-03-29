#ifndef __FORMFACTORSET_HPP
#define __FORMFACTORSET_HPP

/*
 *   Store a collection for FormFactor objects for the nuclei have them for
 */

//#include <FormFactor.hpp>

namespace loos {

    class FormFactorSet
    {

    public:
        FormFactorSet()
            {
            setup();
            }

        FormFactor get(uint i) {
            return _map[i];
        }

    private:
        std::map<uint,FormFactor> _map;
        void setup();
    };
}



#endif
