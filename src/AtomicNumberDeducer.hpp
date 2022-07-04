#if !defined(LOOS_ATOMIC_NUMBER_DEDUCER_HPP)
#define LOOS_ATOMIC_NUMBER_DEDUCER_HPP

#include <vector>
#include <string>


namespace loos {

  namespace internal {

    class AtomicNumberDeducer {
      typedef std::pair<double, unsigned int>   MassNumber;
    public:
      AtomicNumberDeducer() {
        initialize();
      }

      unsigned int deduceFromMass(const double mass, const double tolerance);
      std::string deduceName(const double mass, const double tolerance);

    private:
      void initialize();
      std::vector<MassNumber> element_table;
      std::vector<std::string> name_table;
    };

  };

  //! Deduce an atomic number from the mass
  /**
   * Only the first 96 elements are included in LOOS.  An atomic
   * number of 0 is returned if the mass is not found within LOOS' table.
   */
  unsigned int deduceAtomicNumberFromMass(const double mass, const double tolerance = 0.1);
  std::string deduceElementNameFromMass(const double mass, const double tolerance = 0.1);

};




#endif
