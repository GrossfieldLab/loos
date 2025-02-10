#include <AtomicNumberDeducer.hpp>
#include <cmath>
#include <string>



namespace loos {

  namespace internal {

    unsigned int AtomicNumberDeducer::deduceFromMass(const double mass, const double tolerance) {

      std::vector<MassNumber>::iterator i;
      for (i = element_table.begin(); i != element_table.end(); ++i)
        if (std::abs(i->first - mass) < tolerance)
          break;

      if (i != element_table.end())
        return(i->second);

      return(0);
    }

    double AtomicNumberDeducer::deduceMass(const unsigned int number) {
      static internal::AtomicNumberDeducer deducer;
      std::vector<MassNumber>::iterator i;
      for (i = deducer.element_table.begin(); i != deducer.element_table.end(); ++i)
        if (i->second == number)
          break;
      if (i != element_table.end()) {
        return(i->first);
      }
      return(0);
    }

    std::string AtomicNumberDeducer::deduceName(const double mass,
                                                const double tolerance) {

        unsigned int number = deduceFromMass(mass, tolerance);
        return(name_table[number]);
    }



    void AtomicNumberDeducer::initialize() {

      // These are most frequent based on an analysis of our PSF files
      element_table.push_back(MassNumber(1.008, 1));		// H
      element_table.push_back(MassNumber(15.999, 8));		// O
      element_table.push_back(MassNumber(12.011, 6));		// C
      element_table.push_back(MassNumber(14.007, 7));		// N
      element_table.push_back(MassNumber(30.973762, 15));	// P
      element_table.push_back(MassNumber(35.45, 17));		// Cl
      element_table.push_back(MassNumber(22.98976928, 11));	// Na
      element_table.push_back(MassNumber(32.06, 16));		// S

      // Guesses about what is next most likely
      element_table.push_back(MassNumber(39.0983, 19));		// K
      element_table.push_back(MassNumber(40.078, 20));		// Ca
      element_table.push_back(MassNumber(54.938045, 25));	// Mn
      element_table.push_back(MassNumber(55.845, 26));		// Fe
      element_table.push_back(MassNumber(24.3050, 12));		// Mg
      element_table.push_back(MassNumber(65.38, 30));		// Zn

      // Everything else (hopefully uncommon)
      element_table.push_back(MassNumber(4.002602, 2));		// He
      element_table.push_back(MassNumber(6.94, 3));		// Li
      element_table.push_back(MassNumber(9.012182, 4));		// Be
      element_table.push_back(MassNumber(10.81, 5));		// B
      element_table.push_back(MassNumber(18.9984032, 9));	// F
      element_table.push_back(MassNumber(20.1797, 10));		// Ne
      element_table.push_back(MassNumber(26.9815386, 13));	// Al
      element_table.push_back(MassNumber(28.085, 14));		// Si
      element_table.push_back(MassNumber(39.948, 18));		// Ar
      element_table.push_back(MassNumber(44.955912, 21));	// Sc
      element_table.push_back(MassNumber(47.867, 22));		// Ti
      element_table.push_back(MassNumber(50.9415, 23));		// V
      element_table.push_back(MassNumber(51.9961, 24));		// Cr
      element_table.push_back(MassNumber(58.6934, 28));		// Ni
      element_table.push_back(MassNumber(58.933195, 27));	// Co
      element_table.push_back(MassNumber(63.546, 29));		// Cu
      element_table.push_back(MassNumber(69.723, 31));		// Ga
      element_table.push_back(MassNumber(72.63, 32));		// Ge
      element_table.push_back(MassNumber(74.92160, 33));	// As
      element_table.push_back(MassNumber(78.96, 34));		// Se
      element_table.push_back(MassNumber(79.904, 35));		// Br
      element_table.push_back(MassNumber(83.798, 36));		// Kr
      element_table.push_back(MassNumber(85.4678, 37));		// Rb
      element_table.push_back(MassNumber(87.62, 38));		// Sr
      element_table.push_back(MassNumber(88.90585, 39));	// Y
      element_table.push_back(MassNumber(91.224, 40));		// Zr
      element_table.push_back(MassNumber(92.90638, 41));	// Nb
      element_table.push_back(MassNumber(95.96, 42));		// Mo
      element_table.push_back(MassNumber(98, 43));		// Tc
      element_table.push_back(MassNumber(101.07, 44));		// Ru
      element_table.push_back(MassNumber(102.90550, 45));	// Rh
      element_table.push_back(MassNumber(106.42, 46));		// Pd
      element_table.push_back(MassNumber(107.8682, 47));	// Ag
      element_table.push_back(MassNumber(112.411, 48));		// Cd
      element_table.push_back(MassNumber(114.818, 49));		// In
      element_table.push_back(MassNumber(118.710, 50));		// Sn
      element_table.push_back(MassNumber(121.760, 51));		// Sb
      element_table.push_back(MassNumber(126.90447, 53));	// I
      element_table.push_back(MassNumber(127.60, 52));		// Te
      element_table.push_back(MassNumber(131.293, 54));		// Xe
      element_table.push_back(MassNumber(132.9054519, 55));	// Cs
      element_table.push_back(MassNumber(137.327, 56));		// Ba
      element_table.push_back(MassNumber(138.90547, 57));	// La
      element_table.push_back(MassNumber(140.116, 58));		// Ce
      element_table.push_back(MassNumber(140.90765, 59));	// Pr
      element_table.push_back(MassNumber(144.242, 60));		// Nd
      element_table.push_back(MassNumber(145, 61));		// Pm
      element_table.push_back(MassNumber(150.36, 62));		// Sm
      element_table.push_back(MassNumber(151.964, 63));		// Eu
      element_table.push_back(MassNumber(157.25, 64));		// Gd
      element_table.push_back(MassNumber(158.92535, 65));	// Tb
      element_table.push_back(MassNumber(162.500, 66));		// Dy
      element_table.push_back(MassNumber(164.93032, 67));	// Ho
      element_table.push_back(MassNumber(167.259, 68));		// Er
      element_table.push_back(MassNumber(168.93421, 69));	// Tm
      element_table.push_back(MassNumber(173.054, 70));		// Yb
      element_table.push_back(MassNumber(174.9668, 71));	// Lu
      element_table.push_back(MassNumber(178.49, 72));		// Hf
      element_table.push_back(MassNumber(180.94788, 73));	// Ta
      element_table.push_back(MassNumber(183.84, 74));		// W
      element_table.push_back(MassNumber(186.207, 75));		// Re
      element_table.push_back(MassNumber(190.23, 76));		// Os
      element_table.push_back(MassNumber(192.217, 77));		// Ir
      element_table.push_back(MassNumber(195.084, 78));		// Pt
      element_table.push_back(MassNumber(196.966569, 79));	// Au
      element_table.push_back(MassNumber(200.59, 80));		// Hg
      element_table.push_back(MassNumber(204.38, 81));		// Tl
      element_table.push_back(MassNumber(207.2, 82));		// Pb
      element_table.push_back(MassNumber(208.98040, 83));	// Bi
      element_table.push_back(MassNumber(209, 84));		// Po
      element_table.push_back(MassNumber(210, 85));		// At
      element_table.push_back(MassNumber(222, 86));		// Rn
      element_table.push_back(MassNumber(223, 87));		// Fr
      element_table.push_back(MassNumber(226, 88));		// Ra
      element_table.push_back(MassNumber(227, 89));		// Ac
      element_table.push_back(MassNumber(231.03588, 91));	// Pa
      element_table.push_back(MassNumber(232.03806, 90));	// Th
      element_table.push_back(MassNumber(237, 93));		// Np
      element_table.push_back(MassNumber(238.02891, 92));	// U
      element_table.push_back(MassNumber(243, 95));		// Am
      element_table.push_back(MassNumber(244, 94));		// Pu
      element_table.push_back(MassNumber(247, 96));		// Cm

      name_table.push_back(std::string("NO"));
      name_table.push_back(std::string("H"));
      name_table.push_back(std::string("HE"));
      name_table.push_back(std::string("LI"));
      name_table.push_back(std::string("BE"));
      name_table.push_back(std::string("B"));
      name_table.push_back(std::string("C"));
      name_table.push_back(std::string("N"));
      name_table.push_back(std::string("O"));
      name_table.push_back(std::string("F"));
      name_table.push_back(std::string("NE"));
      name_table.push_back(std::string("NA"));
      name_table.push_back(std::string("MG"));
      name_table.push_back(std::string("AL"));
      name_table.push_back(std::string("SI"));
      name_table.push_back(std::string("P"));
      name_table.push_back(std::string("S"));
      name_table.push_back(std::string("CL"));
      name_table.push_back(std::string("AR"));
      name_table.push_back(std::string("K"));
      name_table.push_back(std::string("CA"));
      name_table.push_back(std::string("SC"));
      name_table.push_back(std::string("TI"));
      name_table.push_back(std::string("V"));
      name_table.push_back(std::string("CR"));
      name_table.push_back(std::string("MN"));
      name_table.push_back(std::string("FE"));
      name_table.push_back(std::string("CO"));
      name_table.push_back(std::string("NI"));
      name_table.push_back(std::string("CU"));
      name_table.push_back(std::string("ZN"));
      name_table.push_back(std::string("GA"));
      name_table.push_back(std::string("GE"));
      name_table.push_back(std::string("AS"));
      name_table.push_back(std::string("SE"));
      name_table.push_back(std::string("BR"));
      name_table.push_back(std::string("KR"));
      name_table.push_back(std::string("RB"));
      name_table.push_back(std::string("SR"));
      name_table.push_back(std::string("Y"));
      name_table.push_back(std::string("ZR"));
      name_table.push_back(std::string("NB"));
      name_table.push_back(std::string("MO"));
      name_table.push_back(std::string("TC"));
      name_table.push_back(std::string("RU"));
      name_table.push_back(std::string("RH"));
      name_table.push_back(std::string("PD"));
      name_table.push_back(std::string("AG"));
      name_table.push_back(std::string("CD"));
      name_table.push_back(std::string("IN"));
      name_table.push_back(std::string("SN"));
      name_table.push_back(std::string("SB"));
      name_table.push_back(std::string("TE"));
      name_table.push_back(std::string("I"));
      name_table.push_back(std::string("XE"));
      name_table.push_back(std::string("CS"));
      name_table.push_back(std::string("BA"));
      name_table.push_back(std::string("LA"));
      name_table.push_back(std::string("CE"));
      name_table.push_back(std::string("PR"));
      name_table.push_back(std::string("ND"));
      name_table.push_back(std::string("PM"));
      name_table.push_back(std::string("SM"));
      name_table.push_back(std::string("EU"));
      name_table.push_back(std::string("GD"));
      name_table.push_back(std::string("TB"));
      name_table.push_back(std::string("DY"));
      name_table.push_back(std::string("HO"));
      name_table.push_back(std::string("ER"));
      name_table.push_back(std::string("TM"));
      name_table.push_back(std::string("YB"));
      name_table.push_back(std::string("LU"));
      name_table.push_back(std::string("HF"));
      name_table.push_back(std::string("TA"));
      name_table.push_back(std::string("W"));
      name_table.push_back(std::string("RE"));
      name_table.push_back(std::string("OS"));
      name_table.push_back(std::string("IR"));
      name_table.push_back(std::string("PT"));
      name_table.push_back(std::string("AU"));
      name_table.push_back(std::string("HG"));
      name_table.push_back(std::string("TL"));
      name_table.push_back(std::string("PB"));
      name_table.push_back(std::string("BI"));
      name_table.push_back(std::string("PO"));
      name_table.push_back(std::string("AT"));
      name_table.push_back(std::string("RN"));
      name_table.push_back(std::string("FR"));
      name_table.push_back(std::string("RA"));
      name_table.push_back(std::string("AC"));
      name_table.push_back(std::string("TH"));
      name_table.push_back(std::string("PA"));
      name_table.push_back(std::string("U"));
      name_table.push_back(std::string("NP"));
      name_table.push_back(std::string("PU"));
      name_table.push_back(std::string("AM"));
      name_table.push_back(std::string("CM"));

    }

  };


  unsigned int deduceAtomicNumberFromMass(const double mass,
                                          const double tolerance) {
    static internal::AtomicNumberDeducer deducer;
    return(deducer.deduceFromMass(mass, tolerance));
  }

  std::string deduceElementNameFromMass(const double mass,
                                        const double tolerance) {
      static internal::AtomicNumberDeducer deducer;
      return deducer.deduceName(mass, tolerance);
  }

  double deduceMassFromAtomicNumber(const unsigned int number) {
    static internal::AtomicNumberDeducer deducer;
    return(deducer.deduceMass(number));
  }
};
