#ifndef __FORMFACTORSET_HPP
#define __FORMFACTORSET_HPP

/*
 *   Store a collection for FormFactor objects for the nuclei have them for
 */

//#include <FormFactor.hpp>

namespace loos {

    class FormFactorSet
    {
    private:
        std::map<uint,FormFactor> _map;
        void setup();

    public:
        FormFactorSet()
            {
            setup();
            }

        FormFactor& operator[](uint i) {
            return _map[i];
        }
    };
}



#endif
