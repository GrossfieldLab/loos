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

#include <loos.hpp>
#include <FormFactor.hpp>

namespace loos {

    // TODO: this should be done with a functor
    double FormFactor::compute(double q) {
        if (atomicNumber <= 7) {
            return smallCompute(q);
        }
        else {
            return bigCompute(q);
        }
    }


    double FormFactor::smallCompute(double q)
        {
        // unpack the coeffs into Szaloki's notation
        // Equation 4
        double a = coeff[0];
        double b1 = coeff[1];
        double c = coeff[2];
        double q1 = coeff[3];
        double b2 = coeff[4];
        double q2 = coeff[5];
        double b3 = coeff[6];
        double q3 = coeff[7];
        double b4 = coeff[8];
        double q4 = coeff[9];

        if ( (q < 0) || (q > q4))
            {
            cerr << "q value out of bounds: ";
            cerr << q << "\t";
            cerr << q4 << "\t";
            exit(-1);
            }

        double val = 0.0;
        if (q <= q1)
            {
            val = f11(q,a,b1,c);
            }
        else if (q <= q2)
            {
            val = f12(q,a,b1,c,b2,q1);
            }
        else if (q <= q3)
            {
            val = f13(q,a,b1,c,b2,q1,b3,q2);
            }
        else if (q <= q4)
            {
            val = f14(q,a,b1,c,b2,q1,b3,q2,b4,q3);
            }

        return val;
        }

    double FormFactor::f11(double q, double a, double b1, double c)
        {
        double val = a*exp(-b1*q) + (atomicNumber - a)*exp(-c*q);
        return val;
        }

    double FormFactor::f12(double q, double a, double b1, double c, double b2, double q1)
        {
        double val = f11(q1,a,b1,c);
        double val2 = exp(b2*(q1 - q));  // sign flipped, Szaloki wrong
        return val*val2;
        }

    double FormFactor::f13(double q, double a, double b1, double c, double b2, double q1, double b3, double q2)
        {
        double val = f12(q2, a, b1, c, b2, q1);
        double val2 = exp(b3*(q2-q));  // sign flipped, Szaloki wrong
        return val*val2;
        }

    double FormFactor::f14(double q, double a, double b1, double c, double b2, double q1, double b3, double q2, double b4, double q3)
        {
        double val = f13(q3,a,b1,c,b2,q1,b3,q2);
        double val2 = pow(q/q3, -b4);  // sign flipped, Szaloki wrong
        return val*val2;
        }

    double FormFactor::bigCompute(double q)
        {
        // unpack the coeffs into Szaloki's notation
        // Equation 5
        double a = coeff[0];
        double b1 = coeff[1];
        double c = coeff[2];
        double q1 = coeff[3];
        double b2 = coeff[4];
        double q2 = coeff[5];
        double b3 = coeff[6];
        double q3 = coeff[7];

        if ( (q < 0) || (q > q3))
            {
            cerr << "q value out of bounds: ";
            cerr << q << "\t";
            cerr << q3 << "\t";
            exit(-1);
            }

        double val = 0.0;
        if (q <= q1)
            {
            val = f21(q,a,b1,c);
            }
        else if (q <= q2)
            {
            val = f22(q,a,b1,c,b2,q1);
            }
        else if (q <= q3)
            {
            val = f23(q,a,b1,c,b2,q1,b3,q2);
            }

        return val;
        }

    double FormFactor::f21(double q, double a, double b1, double c)
        {
        double val = a*exp(-b1*q) + (atomicNumber-a)*exp(-c*q);
        return val;
        }

    double FormFactor::f22(double q, double a, double b1, double c, double b2, double q1)
        {
        double val = f21(q1, a, b1, c);
        double val2 = exp(b2*(q1-q)); // sign flipped, Szaloki wrong
        return val*val2;
        }

    double FormFactor::f23(double q, double a, double b1, double c, double b2, double q1, double b3, double q2) 
        {
        double val = f22(q2, a, b1, c, b2, q1);
        double val2 = exp(b3*(q2-q)); // sign flipped, Szaloki wrong
        return val*val2;
        }


}
