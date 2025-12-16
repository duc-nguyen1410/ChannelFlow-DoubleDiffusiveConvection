/**
 * Solves vector "tau" eqns (Canuto & Hussaini eqn 7.3.18-20)
 * with modification of Dirichlet BC to Neumann BC
 * Original author: Duc Nguyen
 */

#ifndef DDCTAUSOLVER_H
#define DDCTAUSOLVER_H
#include "modules/ddc/macros.h"
#include "cfbasics/cfvector.h"
#include "cfbasics/mathdefs.h"
#include "channelflow/chebyshev.h"
#include "channelflow/helmholtz.h"
#include "modules/ddc/ddchelmholtz.h"
#include "channelflow/tausolver.h"

namespace chflow {
#ifdef FREESLIP
// Class for solving 7.3.18-7.3.20 of Canuto & Hussaini
//      nu u''_jk(y) - lambda u_jk(y) - grad P_jk = -R_jk,
//                                       div u_jk = 0
//                                      u'_jk(+-1) = 0
//
// where u_jk(y) is the vector-valued jkth xz-Fourier coeff of u(x,y,z)
//          P(y) is in R
// and the vector operators are interpreted accordingingly (diff in x
// equal multiplication by 2 pi k).

class DDCTauSolver : public TauSolver {
   public:
   
};
#endif
}  // namespace chflow
#endif  // DDCTAUSOLVER_H