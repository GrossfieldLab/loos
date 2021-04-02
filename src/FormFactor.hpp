#ifndef __FORMFACTOR_HPP__
#define __FORMFACTOR_HPP__
/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


/*
 *  Compute the Fourier space atomic form factors for chemical elements
 *
 *  The new coefficients come from:
 *        Szaloki, X-ray Spectrometry (1996), V25, 21-28
 *        Note: Szaloki screwed up the signs of some of his parameters
 *              (b2, b3, and b4) so my equations don't look exactly like his.
 */

#include <vector>

#define PI 3.1415926535897931

namespace loos {
    class FormFactor
    {
    public:
        std::vector<double> coeff;
        int atomicNumber;

        FormFactor()
            {
            }

        FormFactor(int atomicNumber) : atomicNumber(atomicNumber)
            {

            if (atomicNumber == 1)
                {
                coeff.push_back(3.566);
                coeff.push_back(-1.143);
                coeff.push_back(-2.243);
                coeff.push_back(0.20);
                coeff.push_back(6.102);
                coeff.push_back(0.6);
                coeff.push_back(4.442);
                coeff.push_back(0.90);
                coeff.push_back(3.921);
                coeff.push_back(15.0);
                }
            else if (atomicNumber == 6) // Table 1 from Szaloki, z=6
                {
                coeff.push_back(7.366);
                coeff.push_back(0.745);
                coeff.push_back(-3.209);
                coeff.push_back(0.25);
                coeff.push_back(2.395);
                coeff.push_back(0.5);
                coeff.push_back(1.026);
                coeff.push_back(2.50);
                coeff.push_back(3.258);
                coeff.push_back(8.0);
                }
            else if (atomicNumber == 7) // Table 1 from Szaloki, z=7
                {
                coeff.push_back(8.657);
                coeff.push_back(0.222);
                coeff.push_back(-3.815);
                coeff.push_back(0.25);
                coeff.push_back(2.787);
                coeff.push_back(0.5);
                coeff.push_back(0.878);
                coeff.push_back(2.50);
                coeff.push_back(3.003);
                coeff.push_back(8.0);
                }
            else if (atomicNumber == 8) // Table 2 from Szaloki, z=8
                {
                coeff.push_back(-2.038);
                coeff.push_back(17.634);
                coeff.push_back(2.887);
                coeff.push_back(0.50);
                coeff.push_back(1.339);
                coeff.push_back(0.839);
                coeff.push_back(0.718);
                coeff.push_back(7.0);
                }
            else if (atomicNumber == 15) // Table 2 from Szaloki, z=15
                {
                coeff.push_back(-0.998);
                coeff.push_back(45.579);
                coeff.push_back(2.055);
                coeff.push_back(0.30);
                coeff.push_back(1.425);
                coeff.push_back(1.403);
                coeff.push_back(0.413);
                coeff.push_back(7.0);
                }
            else if (atomicNumber == 16) // Table 2 from Szaloki, z=15
                {
                coeff.push_back(-1.457);
                coeff.push_back(33.964);
                coeff.push_back(2.154);
                coeff.push_back(0.30);
                coeff.push_back(1.321);
                coeff.push_back(1.581);
                coeff.push_back(0.373);
                coeff.push_back(10.0);
                }
            else
                {
                std::cerr << "Unsupported atom type: " << atomicNumber << std::endl;
                exit(-1);
                }
            }

        double compute(double q);
    private:
        double smallCompute(double q);

        double f11(double q, double a, double b1, double c);
        double f12(double q, double a, double b1, double c, double b2, double q1);
        double f13(double q, double a, double b1, double c, double b2, double q1, double b3, double q2);
        double f14(double q, double a, double b1, double c, double b2, double q1, double b3, double q2, double b4, double q3);

        double bigCompute(double q);
        double f21(double q, double a, double b1, double c);
        double f22(double q, double a, double b1, double c, double b2, double q1);
        double f23(double q, double a, double b1, double c, double b2, double q1, double b3, double q2);

    };

}
#endif
