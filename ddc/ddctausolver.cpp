/**
 * Solves vector "tau" eqns (Canuto & Hussaini eqn 7.3.18-20)
 * with modification of Dirichlet BC to Neumann BC
 * Original author: Duc Nguyen
 */

#include "modules/ddc/ddctausolver.h"
#include <iomanip>
#include "cfbasics/mathdefs.h"
namespace chflow {

// Class for solving 7.3.18-7.3.20 of Canuto & Hussaini
//      nu u''_jk(y) - lambda u_jk(y) - grad P_jk = -R_jk,
//                                       div u_jk = 0
//                                      u'_jk(+-1) = 0
//
// where u_jk(y) is the vector-valued jkth xz-Fourier coeff of u(x,y,z)
//          P(y) is in R
// and the vector operators are interpreted accordingingly (diff in x
// equal multiplication by 2 pi k).
#ifdef FREESLIP

int n_func(int k, int Nb) {
    if (k == 0)
        return Nb - 1;
    else if (k == Nb)
        return 0;
    else if (k % 2 == 0)
        return 2 * (Nb - 1);
    else
        return 2 * Nb;
}
#ifndef NDEBUG
const Real MINIMUM_DISCRIMINANT = 1e-4;
#endif
const Real EPSILON = 1e-7;

#endif
}  // namespace chflow
