/**
 * Solution of Helmholtz eqn in Chebyshev expansions
 * EQ.: nu u'' - lambda u = f
 * with modification of Dirichlet BC to Neumann BC
 * Original author: Duc Nguyen
 */


#include "modules/ddc/ddchelmholtz.h"
#include "channelflow/chebyshev.h"
namespace chflow {
#ifdef FREESLIP
const Real EPSILON = 1.0;

void DDCHelmholtzSolver::solve(ChebyCoeff& u, const ChebyCoeff& f, Real ua, Real ub, bool neumann) const {
    assert(f.state() == Spectral);
    Be_.multiplyStrided(f, u, 0, 2);
    Bo_.multiplyStrided(f, u, 1, 2);

    u[0] = (ub + ua) / 2;
    u[1] = (ub - ua) / 2;

    // Solve for A u = g
    Ae_.ULsolveStrided(u, 0, 2);
    Ao_.ULsolveStrided(u, 1, 2);

    //#ifdef DEBUG
    // verify(u, f, ua, ub);
    //#endif
    u.setState(Spectral);
}
#endif
}  // namespace chflow
