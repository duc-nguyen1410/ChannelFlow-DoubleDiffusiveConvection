/**
 * Original author: Duc Nguyen
 */
#include <iostream>
#include <cmath>
using namespace std;
#include "modules/ddc/macros.h"
#include "modules/ddc/boundaryCondition.h"
#include "modules/ddc/dde.h"

namespace chflow {

void momentumNL(const FlowField& u, const FlowField& T, const FlowField& S, 
                ChebyCoeff Ubase, ChebyCoeff Wbase, 
                FlowField& f, FlowField& tmp, DDCFlags flags) {
    // goal: (u*grad)u - P2*(P3*T-P4*S)*(sin(gammax)*ex+cos(gammax)*ey)
    Real Rey = flags.Rey;
    Real Pr = flags.Pr;
    Real Ra = flags.Ra;
    Real Rrho = flags.Rrho;
    Real Rsep = flags.Rsep;
    Real Ri = flags.Ri;
    Real Ek = flags.Ek;
    Real sgammax = sin(flags.gammax);
    Real cgammax = cos(flags.gammax);

    // compute the nonlinear term of NSE in the usual Channelflow style
    navierstokesNL(u, Ubase, Wbase, f, tmp, flags);

    #if defined(P5)||defined(P6)
    // substract the linear temperature+salinity coupling term
    // f -= P2(P3*T-P4*S)*(sin(gammax)*ex+cos(gammax)*ey)
    #ifdef HAVE_MPI
    for (int mz = f.mzlocmin(); mz < f.mzlocmin() + f.Mzloc(); mz++)
        for (int mx = f.mxlocmin(); mx < f.mxlocmin() + f.Mxloc(); mx++)
            for (int ny = 0; ny < f.Ny(); ny++) {
                #if defined(P6) // is for DDC
                f.cmplx(mx, ny, mz, 0) -= P2*(P3*T.cmplx(mx, ny, mz, 0)-P4*S.cmplx(mx, ny, mz, 0))*sgammax;
                f.cmplx(mx, ny, mz, 1) -= P2*(P3*T.cmplx(mx, ny, mz, 0)-P4*S.cmplx(mx, ny, mz, 0))*cgammax;
                #elif defined(P5) // is for pure RBC
                f.cmplx(mx, ny, mz, 0) -= P2*P3*T.cmplx(mx, ny, mz, 0)*sgammax;
                f.cmplx(mx, ny, mz, 1) -= P2*P3*T.cmplx(mx, ny, mz, 0)*cgammax;
                #endif
                #if defined(P8) // P8 is for rotating flows
                // Coriolis contribution: - 1/Ek * sqrt(Pr/Ra) * (ey × u)
                // ey = (0,1,0)  (rotation about vertical y-axis)
                f.cmplx(mx, ny, mz, 0) -= -P8 * u.cmplx(mx, ny, mz, 2); 
                f.cmplx(mx, ny, mz, 1) -= -P8 * (-1) * u.cmplx(mx, ny, mz, 0);
                #endif
            }
    #else
    for (int ny = 0; ny < f.Ny(); ny++)
        for (int mx = f.mxlocmin(); mx < f.mxlocmin() + f.Mxloc(); mx++)
            for (int mz = f.mzlocmin(); mz < f.mzlocmin() + f.Mzloc(); mz++) {
                #if defined(P6)
                f.cmplx(mx, ny, mz, 0) -= P2*(P3*T.cmplx(mx, ny, mz, 0)-P4*S.cmplx(mx, ny, mz, 0))*sgammax;
                f.cmplx(mx, ny, mz, 1) -= P2*(P3*T.cmplx(mx, ny, mz, 0)-P4*S.cmplx(mx, ny, mz, 0))*cgammax;
                #elif defined(P5)
                f.cmplx(mx, ny, mz, 0) -= P2*P3*T.cmplx(mx, ny, mz, 0)*sgammax;
                f.cmplx(mx, ny, mz, 1) -= P2*P3*T.cmplx(mx, ny, mz, 0)*cgammax;
                #endif
                #if defined(P8) // P8 is for rotating flows
                // Coriolis contribution: - 1/Ek * sqrt(Pr/Ra) * (ey × u)
                // ey = (0,1,0)  (rotation about vertical y-axis)
                f.cmplx(mx, ny, mz, 0) -= -P8 * u.cmplx(mx, ny, mz, 2); 
                f.cmplx(mx, ny, mz, 1) -= -P8 * (-1) * u.cmplx(mx, ny, mz, 0);
                #endif
            }
    #endif
    #endif
    
    // dealiasing modes
    if (flags.dealias_xz())
        f.zeroPaddedModes();
}

void temperatureNL(const FlowField& u_, const FlowField& T_, 
                   ChebyCoeff Ubase, ChebyCoeff Wbase, ChebyCoeff Tbase,
                   FlowField& f, FlowField& tmp, DDCFlags flags) {
    FlowField& u = const_cast<FlowField&>(u_);
    FlowField& T = const_cast<FlowField&>(T_);
    // goal: (u*grad)T

    // f += Base;
    for (int ny = 0; ny < u.Ny(); ++ny) {
        if (u.taskid() == u.task_coeff(0, 0)) {
            u.cmplx(0, ny, 0, 0) += Complex(Ubase(ny), 0.0);
            u.cmplx(0, ny, 0, 2) += Complex(Wbase(ny), 0.0);
        }
        if (u.taskid() == u.task_coeff(0, 0))
            T.cmplx(0, ny, 0, 0) += Complex(Tbase(ny), 0.0);
    }
    if (u.taskid() == u.task_coeff(0, 0)) {
        u.cmplx(0, 0, 0, 1) -= Complex(flags.Vsuck, 0.);
    }

    // compute the nonlinearity (temperature advection (u*grad)T ) analogous to the convectiveNL and store in f
    dotgradScalar(u, T, f, tmp);

    // f -= Base;
    for (int ny = 0; ny < u.Ny(); ++ny) {
        if (u.taskid() == u.task_coeff(0, 0)) {
            u.cmplx(0, ny, 0, 0) -= Complex(Ubase(ny), 0.0);
            u.cmplx(0, ny, 0, 2) -= Complex(Wbase(ny), 0.0);
        }
        if (u.taskid() == u.task_coeff(0, 0))
            T.cmplx(0, ny, 0, 0) -= Complex(Tbase(ny), 0.0);
    }
    if (u.taskid() == u.task_coeff(0, 0)) {
        u.cmplx(0, 0, 0, 1) += Complex(flags.Vsuck, 0.);
    }
    

    // dealiasing modes
    if (flags.dealias_xz())
        f.zeroPaddedModes();
}

void salinityNL(const FlowField& u_, const FlowField& T_, const FlowField& S_, 
                ChebyCoeff Ubase, ChebyCoeff Wbase, ChebyCoeff Sbase,
                   FlowField& f, FlowField& tmp, DDCFlags flags) {
    FlowField& u = const_cast<FlowField&>(u_);
    FlowField& T = const_cast<FlowField&>(T_);
    FlowField& S = const_cast<FlowField&>(S_);
    
     // goal: (u*grad)s

    // f += Base;
    for (int ny = 0; ny < u.Ny(); ++ny) {
        if (u.taskid() == u.task_coeff(0, 0)) {
            u.cmplx(0, ny, 0, 0) += Complex(Ubase(ny), 0.0);
            u.cmplx(0, ny, 0, 2) += Complex(Wbase(ny), 0.0);
        }
        if (u.taskid() == u.task_coeff(0, 0))
            S.cmplx(0, ny, 0, 0) += Complex(Sbase(ny), 0.0);
    }
    if (u.taskid() == u.task_coeff(0, 0)) {
        u.cmplx(0, 0, 0, 1) -= Complex(flags.Vsuck, 0.);
    }

    // compute the nonlinearity (temperature advection (u*grad)S ) analogous to the convectiveNL and store in f
    dotgradScalar(u, S, f, tmp);

    #ifdef P7
        for (int mz = f.mzlocmin(); mz < f.mzlocmin() + f.Mzloc(); mz++)
            for (int mx = f.mxlocmin(); mx < f.mxlocmin() + f.Mxloc(); mx++){
                ComplexChebyCoeff Tk_(f.Ny(), f.a(), f.b(), Spectral);
                ComplexChebyCoeff Rtk_(f.Ny(), f.a(), f.b(), Spectral);
                ComplexChebyCoeff Pyk_(f.Ny(), f.a(), f.b(), Spectral);
                // Extract relevant Fourier modes of T: use Pk and Rxk
                for (int ny = 0; ny < f.Ny(); ++ny)
                    Tk_.set(ny, T.cmplx(mx, ny, mz, 0));

                // (1) Put T" into in R. (Pyk_ is used as tmp workspace)
                diff2(Tk_, Rtk_, Pyk_);

                for (int ny = 0; ny < f.Ny(); ny++) {
                    f.cmplx(mx, ny, mz, 0) -= P7*Rtk_[ny];
                }
            }
    #endif

    // f -= Base;
    for (int ny = 0; ny < u.Ny(); ++ny) {
        if (u.taskid() == u.task_coeff(0, 0)) {
            u.cmplx(0, ny, 0, 0) -= Complex(Ubase(ny), 0.0);
            u.cmplx(0, ny, 0, 2) -= Complex(Wbase(ny), 0.0);
        }
        if (u.taskid() == u.task_coeff(0, 0))
            S.cmplx(0, ny, 0, 0) -= Complex(Sbase(ny), 0.0);
    }
    if (u.taskid() == u.task_coeff(0, 0)) {
        u.cmplx(0, 0, 0, 1) += Complex(flags.Vsuck, 0.);
    }

    // dealiasing modes
    if (flags.dealias_xz())
        f.zeroPaddedModes();
}

DDE::DDE(const std::vector<FlowField>& fields, const DDCFlags& flags)
    : NSE(fields, flags),
      heatsolver_(0),  // heatsolvers are allocated when reset_lambda is called for the first time
      saltsolver_(0),
      flags_(flags),
      Tbase_(),
      Tbaseyy_(),
      Sbase_(),
      Sbaseyy_(),
      Pbasey_(),
      Cu_(),
      Cw_(),
      Ct_(),
      Cs_(),
      nonzCu_(false),
      nonzCw_(false),
      nonzCt_(false),
      nonzCs_(false),
      Tk_(Nyd_, a_, b_, Spectral),
      Rtk_(Nyd_, a_, b_, Spectral),
      Sk_(Nyd_, a_, b_, Spectral),
      Rsk_(Nyd_, a_, b_, Spectral),
      baseflow_(false),
      constraint_(false) {
    
    // set member variables for base flow
    createDDCBaseFlow();
    Pbasey_ = hydrostaticPressureGradientY(Tbase_, Sbase_, flags_);

    // set member variables for contraint
    initDDCConstraint(fields[0]);

    // define the constant terms of the DDE
    createConstants();
}

DDE::DDE(const std::vector<FlowField>& fields, const std::vector<ChebyCoeff>& base, const DDCFlags& flags)
    : NSE(fields, base, flags),
      heatsolver_(0),  // heatsolvers are allocated when reset_lambda is called for the first time
      saltsolver_(0),
      flags_(flags),
      Tbase_(base[2]),
      Tbaseyy_(),
      Sbase_(base[3]),
      Sbaseyy_(),
      Pbasey_(),
      Cu_(),
      Cw_(),
      Ct_(),
      Cs_(),
      nonzCu_(false),
      nonzCw_(false),
      nonzCt_(false),
      nonzCs_(false),
      Tk_(Nyd_, a_, b_, Spectral),
      Rtk_(Nyd_, a_, b_, Spectral),
      Sk_(Nyd_, a_, b_, Spectral),
      Rsk_(Nyd_, a_, b_, Spectral),
      baseflow_(false),
      constraint_(false) {
    

    baseflow_ = true;  // base flow is passed to constructor
    Pbasey_ = hydrostaticPressureGradientY(Tbase_, Sbase_, flags_);

    // set member variables for contraint
    initDDCConstraint(fields[0]);

    // define the constant terms of the DDE
    createConstants();
}

DDE::~DDE() {
    if (tausolver_) {
        for (uint j = 0; j < lambda_t_.size(); ++j) {
            for (int mx = 0; mx < Mxloc_; ++mx) {
                delete[] tausolver_[j][mx];  // undo new #3
                tausolver_[j][mx] = 0;
            }
            delete[] tausolver_[j];  // undo new #2
            tausolver_[j] = 0;
        }
        delete[] tausolver_;  // undo new #1
        tausolver_ = 0;
    }
    #ifdef P5
    if (heatsolver_) {
        for (uint j = 0; j < lambda_t_.size(); ++j) {
            for (int mx = 0; mx < Mxloc_; ++mx) {
                delete[] heatsolver_[j][mx];  // undo new #3
                heatsolver_[j][mx] = 0;
            }
            delete[] heatsolver_[j];  // undo new #2
            heatsolver_[j] = 0;
        }
        delete[] heatsolver_;  // undo new #1
        heatsolver_ = 0;
    }
    #endif
    #ifdef P6
    if (saltsolver_) {
        for (uint j = 0; j < lambda_t_.size(); ++j) {
            for (int mx = 0; mx < Mxloc_; ++mx) {
                delete[] saltsolver_[j][mx];  // undo new #3
                saltsolver_[j][mx] = 0;
            }
            delete[] saltsolver_[j];  // undo new #2
            saltsolver_[j] = 0;
        }
        delete[] saltsolver_;  // undo new #1
        saltsolver_ = 0;
    }
    #endif
}

void DDE::nonlinear(const std::vector<FlowField>& infields, std::vector<FlowField>& outfields) {//infields=[u,T,S,p] and outfields[u,T,S]
    // The first entry in vector must be velocity FlowField, the second a temperature FlowField, and third is salinity FlowField.
    // Pressure as third entry in in/outfields is not touched.
    momentumNL(infields[0], infields[1], infields[2], Ubase_,Wbase_, outfields[0], tmp_, flags_);
    #ifdef P5
    temperatureNL(infields[0], infields[1], Ubase_,Wbase_,Tbase_, outfields[1], tmp_, flags_);
    #endif
    #ifdef P6
    salinityNL(infields[0], infields[1], infields[2], Ubase_,Wbase_,Sbase_, outfields[2], tmp_, flags_);
    #endif
}

void DDE::linear(const std::vector<FlowField>& infields, std::vector<FlowField>& outfields) {//infields=[u,T,S,p] and outfields[u,T,S]
    // Method takes input fields {u,T,S,press} and computes the linear terms for velocity output field {u,T,S}

    // Make sure user does not expect a pressure output. Outfields should be created outside DDE with DDE::createRHS()
    assert(infields.size() == (outfields.size() + 1));
    const int kxmax = infields[0].kxmax();
    const int kzmax = infields[0].kzmax();
    const Real Rey = flags_.Rey;
    const Real Pr = flags_.Pr;
    const Real Ra = flags_.Ra;
    const Real Le = flags_.Le;
    const Real Rrho = flags_.Rrho;
    const Real Ri = flags_.Ri;
    

    // Loop over Fourier modes. 2nd derivative and summation of linear term
    // is most sufficient on ComplexChebyCoeff. Therefore, the old loop structure is kept.
    for (lint mx = mxlocmin_; mx < mxlocmin_ + Mxloc_; ++mx) {
        const int kx = infields[0].kx(mx);
        for (lint mz = mzlocmin_; mz < mzlocmin_ + Mzloc_; ++mz) {
            const int kz = infields[0].kz(mz);

            // Skip last and aliased modes
            if ((kx == kxmax || kz == kzmax) || (flags_.dealias_xz() && isAliasedMode(kx, kz)))
                break;

            // Compute linear momentum terms
            //================================
            // Goal is to compute
            // L = Pr*u" - Pr*kappa^2*u - grad(p) [+ d2y2Ubase - grad(p0)]

            // Extract relevant Fourier modes of uj and qj
            // get (u) into u_(uk_,vk_,wk_) and get p_
            for (int ny = 0; ny < Nyd_; ++ny) {
                uk_.set(ny, infields[0].cmplx(mx, ny, mz, 0)); // uk = u
                vk_.set(ny, infields[0].cmplx(mx, ny, mz, 1)); // vk = v
                wk_.set(ny, infields[0].cmplx(mx, ny, mz, 2)); // wk = w
                Pk_.set(ny, infields[3].cmplx(mx, ny, mz, 0));  // pressure is 4. entry in fields, Pk = p
            }

            // (1) Put nu uj" into in R. (Pyk_ is used as tmp workspace)
            // get (\nu*u"=(\nu*u)") into variable R(Ruk_,Rvk_,Rwk_)
            // diff2(in,out,temp) --> out = d2(in)/dy2
            diff2(uk_, Ruk_, Pyk_); // Ruk = d2(u)/dy2
            diff2(vk_, Rvk_, Pyk_); // Rvk = d2(v)/dy2
            diff2(wk_, Rwk_, Pyk_); // Rwk = d2(w)/dy2

            // (2) Put qn' into Pyk (compute y-comp of pressure gradient).
            // get grad(p(y))=(dx*p_=Dx*pk_,dy*p_=dpk_/dy,dz*p_=Dz*pk_)
            // diff(in,out) --> out=d(in)/dy
            diff(Pk_, Pyk_); // Pk = dp/dy

            // (3) Summation of all derivative terms and assignment to output velocity field.
            const Real k2 = 4 * pi * pi * (square(kx / Lx_) + square(kz / Lz_));
            const Complex Dx = infields[0].Dx(mx); // get dx
            const Complex Dz = infields[0].Dz(mz); // get dz
            for (int ny = 0; ny < Nyd_; ++ny) {
                outfields[0].cmplx(mx, ny, mz, 0) = P1*Ruk_[ny] - P1*k2*uk_[ny] - Dx * Pk_[ny]; // Pr*d2(u)/dy2-Pr*nablaxz^2*(u)-dp/dx
                outfields[0].cmplx(mx, ny, mz, 1) = P1*Rvk_[ny] - P1*k2*vk_[ny] - Pyk_[ny];     // Pr*d2(v)/dy2-Pr*nablaxz^2*(v)-dp/dy
                outfields[0].cmplx(mx, ny, mz, 2) = P1*Rwk_[ny] - P1*k2*wk_[ny] - Dz * Pk_[ny]; // Pr*d2(w)/dy2-Pr*nablaxz^2*(w)-dp/dz
            }
            // (4) Add const. terms Cu=d2y2Ubase-dPdx
            if (kx == 0 && kz == 0) {
                // L includes const dissipation term of Ubase and Wbase: nu Uyy, nu Wyy
                if (Ubaseyy_.length() > 0)
                    for (int ny = 0; ny < My_; ++ny)
                        outfields[0].cmplx(mx, ny, mz, 0) += Complex(P1 * Ubaseyy_[ny], 0);
                if (Wbaseyy_.length() > 0)
                    for (int ny = 0; ny < My_; ++ny)
                        outfields[0].cmplx(mx, ny, mz, 2) += Complex(P1 * Wbaseyy_[ny], 0);

                // Add base pressure gradient depending on the constraint
                if (flags_.constraint == PressureGradient) {
                    // dPdx [is grad(p0)] is supplied as dPdxRef []
                    outfields[0].cmplx(mx, 0, mz, 0) -= Complex(dPdxRef_, 0);
                    outfields[0].cmplx(mx, 0, mz, 2) -= Complex(dPdzRef_, 0);
                } else {  // const bulk velocity [is pressure gradient grad(p0) constructed from input velocity field]
                    // actual dPdx is unknown but defined by constraint of bulk velocity
                    // Determine actual dPdx from Ubase + u.
                    Real Ly = b_ - a_;
                    diff(uk_, Ruk_);
                    diff(wk_, Rwk_);
                    Real dPdxAct = Re(Ruk_.eval_b() - Ruk_.eval_a()) / Ly;
                    Real dPdzAct = Re(Rwk_.eval_b() - Rwk_.eval_a()) / Ly;
                    ChebyCoeff Ubasey = diff(Ubase_);
                    ChebyCoeff Wbasey = diff(Wbase_);
                    if (Ubase_.length() != 0)
                        dPdxAct += P1 * (Ubasey.eval_b() - Ubasey.eval_a()) / Ly;
                    if (Wbase_.length() != 0)
                        dPdzAct += P1 * (Wbasey.eval_b() - Wbasey.eval_a()) / Ly;
                    // add press. gradient to linear term
                    outfields[0].cmplx(mx, 0, mz, 0) -= Complex(dPdxAct, 0);
                    outfields[0].cmplx(mx, 0, mz, 2) -= Complex(dPdzAct, 0);
                }
            }  // End of const. terms

            #ifdef P5
            // Compute linear heat equation terms
            //================================
            // Goal is to compute
            // L = T" - kappa^2*u [+ d2y2Tbase]

            // Extract relevant Fourier modes of T: use Pk and Rxk
            for (int ny = 0; ny < Nyd_; ++ny)
                Tk_.set(ny, infields[1].cmplx(mx, ny, mz, 0));

            // (1) Put kappa T" into in R. (Pyk_ is used as tmp workspace)
            diff2(Tk_, Rtk_, Pyk_);

            // (2+3) Summation of both diffusive terms and linear advection of temperature.
            // k2 and Dx were defined above for the current Fourier mode
            for (int ny = 0; ny < Nyd_; ++ny)
                // from nonlinear:  - Dx*Tk_[ny]*Complex(Ubase_[ny],0)/kappa;
                outfields[1].cmplx(mx, ny, mz, 0) = P5*Rtk_[ny] - P5*k2 * Tk_[ny];

            // (4) Add const. terms [don't need this, because Ct=d2y2Tbase=0]
            if (kx == 0 && kz == 0) {
                // L includes const diffusion term of Tbase: Tbase_yy
                if (nonzCt_)
                    for (int ny = 0; ny < My_; ++ny)
                        outfields[1].cmplx(mx, ny, mz, 0) += Ct_[ny];
            }
            #endif

            #ifdef P6
            // Compute linear salt equation terms
            //================================
            // Goal is to compute
            // L = 1/Le*S" - kappa^2 *1/Le*S [+ d2y2Sbase]

            // Extract relevant Fourier modes of S: use Sk and Rsk
            for (int ny = 0; ny < Nyd_; ++ny)
                Sk_.set(ny, infields[2].cmplx(mx, ny, mz, 0));

            // (1) Put S" into in R. (Pyk_ is used as tmp workspace)
            diff2(Sk_, Rsk_, Pyk_);

            // (2+3) Summation of both diffusive terms and linear advection of temperature.
            // k2 and Dx were defined above for the current Fourier mode
            for (int ny = 0; ny < Nyd_; ++ny)
                // from nonlinear:  - Dx*Tk_[ny]*Complex(Ubase_[ny],0)/kappa;
                outfields[2].cmplx(mx, ny, mz, 0) =  P6*Rsk_[ny] -  P6*k2*Sk_[ny];

            // (4) Add const. terms [don't need this, because Cs=d2y2Sbase=0]
            if (kx == 0 && kz == 0) {
                // L includes const diffusion term of Sbase:  1/Le * Sbase_yy
                if (nonzCs_)
                    for (int ny = 0; ny < My_; ++ny)
                        outfields[2].cmplx(mx, ny, mz, 0) += Cs_[ny];
            }
            #endif
        }
    }  // End of loop over Fourier modes
}

void DDE::solve(std::vector<FlowField>& outfields, const std::vector<FlowField>& rhs, const int s) {
    // Method takes a right hand side {u} and solves for output fields {u,press}
    
    // Make sure user provides correct RHS which can be created outside NSE with NSE::createRHS()
    assert(outfields.size() == (rhs.size() + 1));
    const int kxmax = outfields[0].kxmax();
    const int kzmax = outfields[0].kzmax();

    // Update each Fourier mode with solution of the implicit problem
    for (lint ix = 0; ix < Mxloc_; ++ix) {
        const lint mx = ix + mxlocmin_;
        const int kx = outfields[0].kx(mx);

        for (lint iz = 0; iz < Mzloc_; ++iz) {
            const lint mz = iz + mzlocmin_;
            const int kz = outfields[0].kz(mz);

            // Skip last and aliased modes
            if ((kx == kxmax || kz == kzmax) || (flags_.dealias_xz() && isAliasedMode(kx, kz)))
                break;

            // Construct ComplexChebyCoeff
            for (int ny = 0; ny < Nyd_; ++ny) {
                Ruk_.set(ny, rhs[0].cmplx(mx, ny, mz, 0));
                Rvk_.set(ny, rhs[0].cmplx(mx, ny, mz, 1));
                Rwk_.set(ny, rhs[0].cmplx(mx, ny, mz, 2));
                // negative RHS because HelmholtzSolver solves the negative problem
                #ifdef P5
                Rtk_.set(ny, -rhs[1].cmplx(mx, ny, mz, 0));
                #endif
                #ifdef P6
                Rsk_.set(ny, -rhs[2].cmplx(mx, ny, mz, 0));
                #endif
            }

            
            // Solve the tau equations for momentum
            //=============================
            if (kx != 0 || kz != 0)
                tausolver_[s][ix][iz].solve(uk_, vk_, wk_, Pk_, Ruk_, Rvk_, Rwk_);
            // 	solve(ix,iz,uk_,vk_,wk_,Pk_, Ruk_,Rvk_,Rwk_);
            else {  // kx,kz == 0,0
                // LHS includes also the constant terms C which can be added to RHS
                if (nonzCu_ || nonzCw_) {
                    for (int ny = 0; ny < My_; ++ny) {
                        Ruk_.re[ny] += Cu_.re[ny];
                        Rwk_.re[ny] += Cw_.re[ny];
                    }
                }
                
                
                if (flags_.constraint == PressureGradient) {
                    // pressure is supplied, put on RHS of tau eqn
                    Ruk_.re[0] -= dPdxRef_;
                    Rwk_.re[0] -= dPdzRef_;
                    tausolver_[s][ix][iz].solve(uk_, vk_, wk_, Pk_, Ruk_, Rvk_, Rwk_);
                    // 	  solve(ix,iz,uk_, vk_, wk_, Pk_, Ruk_,Rvk_,Rwk_);
                    // Bulk vel is free variable determined from soln of tau eqn //TODO: write method that computes
                    // UbulkAct everytime it is needed

                } else {  // const bulk velocity
                    // bulk velocity is supplied, use alternative tau solver

                    // Use tausolver with additional variable and constraint:
                    // free variable: dPdxAct at next time-step,
                    // constraint:    UbulkBase + mean(u) = UbulkRef.
                    tausolver_[s][ix][iz].solve(uk_, vk_, wk_, Pk_, dPdxAct_, dPdzAct_, Ruk_, Rvk_, Rwk_,
                                                UbulkRef_ - UbulkBase_, WbulkRef_ - WbulkBase_);
                    // 	  solve(ix,iz,uk_, vk_, wk_, Pk_, dPdxAct_, dPdzAct_,
                    // 				    Ruk_, Rvk_, Rwk_,
                    // 				    UbulkRef_ - UbulkBase_,
                    // 				    WbulkRef_ - WbulkBase_);

                    // test if UbulkRef == UbulkAct = UbulkBase_ + uk_.re.mean()
                    assert((UbulkRef_ - UbulkBase_ - uk_.re.mean()) < 1e-15);
                    // test if WbulkRef == WbulkAct = WbulkBase_ + wk_.re.mean()
                    assert((WbulkRef_ - WbulkBase_ - wk_.re.mean()) < 1e-15);
                }
                    
            }
            // Load solutions into u and p.
            // Because of FFTW complex symmetries
            // The 0,0 mode must be real.
            // For Nx even, the kxmax,0 mode must be real
            // For Nz even, the 0,kzmax mode must be real
            // For Nx,Nz even, the kxmax,kzmax mode must be real
            if ((kx == 0 && kz == 0) || (outfields[0].Nx() % 2 == 0 && kx == kxmax && kz == 0) ||
                (outfields[0].Nz() % 2 == 0 && kz == kzmax && kx == 0) ||
                (outfields[0].Nx() % 2 == 0 && outfields[0].Nz() % 2 == 0 && kx == kxmax && kz == kzmax)) {
                for (int ny = 0; ny < Nyd_; ++ny) {
                    outfields[0].cmplx(mx, ny, mz, 0) = Complex(Re(uk_[ny]), 0.0);
                    outfields[0].cmplx(mx, ny, mz, 1) = Complex(Re(vk_[ny]), 0.0);
                    outfields[0].cmplx(mx, ny, mz, 2) = Complex(Re(wk_[ny]), 0.0);
                    outfields[3].cmplx(mx, ny, mz, 0) = Complex(Re(Pk_[ny]), 0.0);
                }
            }
            // The normal case, for general kx,kz
            else
                for (int ny = 0; ny < Nyd_; ++ny) {
                    outfields[0].cmplx(mx, ny, mz, 0) = uk_[ny];
                    outfields[0].cmplx(mx, ny, mz, 1) = vk_[ny];
                    outfields[0].cmplx(mx, ny, mz, 2) = wk_[ny];
                    outfields[3].cmplx(mx, ny, mz, 0) = Pk_[ny];
                }

            #ifdef P5
            // Solve the helmholtz problem for the heat equation
            //=============================
            if (kx != 0 || kz != 0) {
                heatsolver_[s][ix][iz].solve(Tk_.re, Rtk_.re, 0, 0);  // BC are considered through the base profile
                heatsolver_[s][ix][iz].solve(Tk_.im, Rtk_.im, 0, 0);
            } else {  // kx,kz == 0,0
                // LHS includes also the constant term C=kappa Tbase_yy, which can be added to RHS
                if (nonzCt_) {
                    for (int ny = 0; ny < My_; ++ny)
                        Rtk_.re[ny] -= Ct_.re[ny];
                }

                heatsolver_[s][ix][iz].solve(Tk_.re, Rtk_.re, 0, 0);
                heatsolver_[s][ix][iz].solve(Tk_.im, Rtk_.im, 0, 0);
            }

            // Load solution into T.
            // Because of FFTW complex symmetries
            // The 0,0 mode must be real.
            // For Nx even, the kxmax,0 mode must be real
            // For Nz even, the 0,kzmax mode must be real
            // For Nx,Nz even, the kxmax,kzmax mode must be real
            if ((kx == 0 && kz == 0) || (outfields[1].Nx() % 2 == 0 && kx == kxmax && kz == 0) ||
                (outfields[1].Nz() % 2 == 0 && kz == kzmax && kx == 0) ||
                (outfields[1].Nx() % 2 == 0 && outfields[1].Nz() % 2 == 0 && kx == kxmax && kz == kzmax)) {
                for (int ny = 0; ny < Nyd_; ++ny)
                    outfields[1].cmplx(mx, ny, mz, 0) = Complex(Re(Tk_[ny]), 0.0);

            }
            // The normal case, for general kx,kz
            else
                for (int ny = 0; ny < Nyd_; ++ny)
                    outfields[1].cmplx(mx, ny, mz, 0) = Tk_[ny];
            #endif
            
            #ifdef P6
            // Solve the helmholtz problem for the salt equation
            //=============================
            if (kx != 0 || kz != 0) {
                saltsolver_[s][ix][iz].solve(Sk_.re, Rsk_.re, 0, 0);  // BC are considered through the base profile
                saltsolver_[s][ix][iz].solve(Sk_.im, Rsk_.im, 0, 0);
            } else {  // kx,kz == 0,0
                // LHS includes also the constant term C=1/Le Sbase_yy, which can be added to RHS
                if (nonzCs_) {
                    for (int ny = 0; ny < My_; ++ny)
                        Rsk_.re[ny] -= Cs_.re[ny];
                }

                saltsolver_[s][ix][iz].solve(Sk_.re, Rsk_.re, 0, 0);
                saltsolver_[s][ix][iz].solve(Sk_.im, Rsk_.im, 0, 0);
            }

            // Load solution into S.
            // Because of FFTW complex symmetries
            // The 0,0 mode must be real.
            // For Nx even, the kxmax,0 mode must be real
            // For Nz even, the 0,kzmax mode must be real
            // For Nx,Nz even, the kxmax,kzmax mode must be real
            if ((kx == 0 && kz == 0) || (outfields[2].Nx() % 2 == 0 && kx == kxmax && kz == 0) ||
                (outfields[2].Nz() % 2 == 0 && kz == kzmax && kx == 0) ||
                (outfields[2].Nx() % 2 == 0 && outfields[2].Nz() % 2 == 0 && kx == kxmax && kz == kzmax)) {
                for (int ny = 0; ny < Nyd_; ++ny)
                    outfields[2].cmplx(mx, ny, mz, 0) = Complex(Re(Sk_[ny]), 0.0);

            }
            // The normal case, for general kx,kz
            else
                for (int ny = 0; ny < Nyd_; ++ny)
                    outfields[2].cmplx(mx, ny, mz, 0) = Sk_[ny];
            #endif

        }  // End of loop over Fourier modes
    }
    
}

void DDE::reset_lambda(std::vector<Real> lambda_t) {
    lambda_t_ = lambda_t;
    if (tausolver_ == 0) {  // TauSolver need to be constructed
        // Allocate memory for [Nsubsteps x Mx_ x Mz_] Tausolver cfarray
        tausolver_ = new TauSolver**[lambda_t.size()];  // new #1
        for (uint j = 0; j < lambda_t.size(); ++j) {
            tausolver_[j] = new TauSolver*[Mxloc_];  // new #2
            for (int mx = 0; mx < Mxloc_; ++mx)
                tausolver_[j][mx] = new TauSolver[Mzloc_];  // new #3
        }
    }
    #ifdef P5
    if (heatsolver_ == 0) {
        heatsolver_ = new HelmholtzSolver**[lambda_t.size()];  // new #1
        for (uint j = 0; j < lambda_t.size(); ++j) {
            heatsolver_[j] = new HelmholtzSolver*[Mxloc_];  // new #2
            for (int mx = 0; mx < Mxloc_; ++mx)
                heatsolver_[j][mx] = new HelmholtzSolver[Mzloc_];  // new #3
        }
    }
    #endif
    #ifdef P6
    if (saltsolver_ == 0) {
        saltsolver_ = new HelmholtzSolver**[lambda_t.size()];  // new #1
        for (uint j = 0; j < lambda_t.size(); ++j) {
            saltsolver_[j] = new HelmholtzSolver*[Mxloc_];  // new #2
            for (int mx = 0; mx < Mxloc_; ++mx)
                saltsolver_[j][mx] = new HelmholtzSolver[Mzloc_];  // new #3
        }
    }
    #endif

    // Configure tausolvers
    //   FlowField u=fields_[0];

    const Real Rey = flags_.Rey;
    const Real Pr = flags_.Pr;
    const Real Le = flags_.Le;
    const Real Ra = flags_.Ra;
    const Real Rrho = flags_.Rrho;
    const Real Ri = flags_.Ri;
    // const Real Tau = flags_.Tau;
    const Real c = 4.0 * square(pi);
    //   const int kxmax = u.kxmax();
    //   const int kzmax = u.kzmax();
    for (uint j = 0; j < lambda_t.size(); ++j) {
        for (int mx = 0; mx < Mxloc_; ++mx) {
            int kx = kxloc_[mx];
            for (int mz = 0; mz < Mzloc_; ++mz) {
                int kz = kzloc_[mz];
                Real lambda_tau = lambda_t[j] + P1 * c * (square(kx / Lx_) + square(kz / Lz_));
                #ifdef P5
                Real lambda_heat = lambda_t[j] + P5 * c * (square(kx / Lx_) + square(kz / Lz_));
                #endif
                #ifdef P6
                Real lambda_salt = lambda_t[j] + P6 * c * (square(kx / Lx_) + square(kz / Lz_));
                #endif
                if ((kx != kxmax_ || kz != kzmax_) && (!flags_.dealias_xz() || !isAliasedMode(kx, kz))) {
                    tausolver_[j][mx][mz] =
                        TauSolver(kx, kz, Lx_, Lz_, a_, b_, lambda_tau, P1, Nyd_, flags_.taucorrection);
                    #ifdef P5
                    heatsolver_[j][mx][mz] = HelmholtzSolver(Nyd_, a_, b_, lambda_heat, P5);
                    #endif
                    #ifdef P6
                    saltsolver_[j][mx][mz] = HelmholtzSolver(Nyd_, a_, b_, lambda_salt, P6);
                    #endif
                }
            }
        }
    }
}


const ChebyCoeff& DDE::Ubase() const { return Ubase_; }
const ChebyCoeff& DDE::Wbase() const { return Wbase_; }
const ChebyCoeff& DDE::Tbase() const { return Tbase_; }
const ChebyCoeff& DDE::Sbase() const { return Sbase_; }

std::vector<FlowField> DDE::createRHS(const std::vector<FlowField>& fields) const {
    return {fields[0], fields[1], fields[2]};  // return only velocity, temperature, and salinity fields
}

std::vector<cfarray<FieldSymmetry>> DDE::createSymmVec() const {
    // velocity symmetry
    cfarray<FieldSymmetry> usym = SymmetryList(1);  // unity operation
    if (flags_.symmetries.length() > 0)
        usym = flags_.symmetries;
    // temperature symmetry
    cfarray<FieldSymmetry> tsym = SymmetryList(1);  // unity operation
    if (flags_.tempsymmetries.length() > 0)
        tsym = flags_.tempsymmetries;
    // salinity symmetry
    cfarray<FieldSymmetry> ssym = SymmetryList(1);  // unity operation
    if (flags_.saltsymmetries.length() > 0)
        ssym = flags_.saltsymmetries;
    // pressure symmetry
    cfarray<FieldSymmetry> psym = SymmetryList(1);  // unity operation
    return {usym, tsym, ssym, psym};
}

void DDE::createDDCBaseFlow() {
    Real ulowerwall = flags_.ulowerwall;
    Real uupperwall = flags_.uupperwall;
    Real wlowerwall = flags_.wlowerwall;
    Real wupperwall = flags_.wupperwall;

    switch (flags_.baseflow) {
        case ZeroBase:
            Ubase_ = ChebyCoeff(My_, a_, b_, Spectral);
            Wbase_ = ChebyCoeff(My_, a_, b_, Spectral);
            Tbase_ = ChebyCoeff(My_, a_, b_, Spectral);
            Sbase_ = ChebyCoeff(My_, a_, b_, Spectral);
            break;
        case LinearBase:
            std::cerr << "error in DDE::createBaseFlow :\n";
            std::cerr << "LinearBase is not defined in DDC.\n";
            break;
        case ParabolicBase:
            std::cerr << "error in DDE::createBaseFlow :\n";
            std::cerr << "ParabolicBase is not defined in DDC.\n";
            break;
        case SuctionBase:
            std::cerr << "error in DDE::createBaseFlow :\n";
            std::cerr << "ParabolicBase is not defined in DDC.\n";
            break;
        case LaminarBase:
            Ubase_ = laminarVelocityProfile(flags_.gammax, flags_.dPdx, flags_.Ubulk, ulowerwall, uupperwall, a_, b_, My_, flags_);
            Wbase_ = laminarVelocityProfile(0.0, flags_.dPdz, flags_.Wbulk, wlowerwall, wupperwall, a_, b_, My_, flags_);

            Tbase_ = linearTemperatureProfile(a_, b_, My_, flags_);
            Sbase_ = linearSalinityProfile(a_, b_, My_, flags_);

            break;
        case ArbitraryBase:
            std::cerr << "error in NSE::createBaseFlow :\n";
            std::cerr << "flags.baseflow is ArbitraryBase.\n";
            std::cerr << "Please provide {Ubase, Wbase} when constructing DNS.\n";
            cferror("");
        default:
            std::cerr << "error in NSE::createBaseFlow :\n";
            std::cerr << "flags.baseflow should be ZeroBase, LinearBase, ParabolicBase, LaminarBase, SuctionBase.\n";
            std::cerr << "Other cases require use of the DNS::DNS(fields, base, flags) constructor.\n";
            cferror("");
    }
    baseflow_ = true;
}

void DDE::createConstants() {
    if ((!baseflow_) || (!constraint_))
        std::cerr << "DDE::createConstants: Not all base flow members have been initialized." << std::endl;

    ComplexChebyCoeff c(My_, a_, b_, Spectral);
    Real Rey = flags_.Rey;
    Real Pr = flags_.Pr;
    Real Le = flags_.Le;
    Real Ra = flags_.Ra;
    Real Rrho = flags_.Rrho;
    Real Ri = flags_.Ri;

    // constant u-term
    if (Ubaseyy_.length() > 0) {
        c[0] = Complex(P1 * Ubaseyy_[0], 0);
        if ((abs(c[0]) > 1e-15) && !nonzCu_)
            nonzCu_ = true;
        for (int ny = 1; ny < My_; ++ny) {
            c[ny] = Complex(P1 * Ubaseyy_[ny], 0);
            if ((abs(c[ny]) > 1e-15) && !nonzCu_)
                nonzCu_ = true;
        }
    }
    if (nonzCu_) {
        Cu_ = c;
        *flags_.logstream << "DDC with nonzero U-const." << std::endl;
    }

    // check that constant v-term is zero:
    // Real hydrostaty = - Pbasey_[0];
    // if (abs(hydrostaty) > 1e-15)
    //     std::cerr << "Wall-normal hydrostatic pressure is unballanced" << std::endl;
    // for (int ny = 1; ny < My_; ++ny)
    //     if (abs(- Pbasey_[ny]) > 1e-15)
    //         std::cerr << "Wall-normal hydrostatic pressure is unballanced" << std::endl;

    // constant w-term:
    if (Wbaseyy_.length() > 0) {
        for (int ny = 0; ny < My_; ++ny) {
            c[ny] = Complex(P1 * Wbaseyy_[ny], 0);
            if ((abs(c[ny]) > 1e-15) && !nonzCw_)
                nonzCw_ = true;
        }
    }
    if (nonzCw_) {
        Cw_ = c;
        *flags_.logstream << "DDC with nonzero W-const." << std::endl;
    }

    #ifdef P5
    // constant t-term:
    if (Tbaseyy_.length() > 0) {
        for (int ny = 1; ny < My_; ++ny) {
            c[ny] = Complex(P5*Tbaseyy_[ny], 0);
            if ((abs(c[ny]) > 1e-15) && !nonzCt_)
                nonzCt_ = true;
        }
    }
    if (nonzCt_) {
        Ct_ = c;
        *flags_.logstream << "DDC with nonzero T-const." << std::endl;
    }
    #endif

    #ifdef P6
    // constant s-term:
    if (Sbaseyy_.length() > 0) {
        for (int ny = 1; ny < My_; ++ny) {
            c[ny] = Complex(P6*Sbaseyy_[ny]
            #ifdef P7
                + P7*Tbaseyy_[ny]
            #endif
            , 0);
            if ((abs(c[ny]) > 1e-15) && !nonzCs_)
                nonzCs_ = true;
        }
    }
    if (nonzCs_) {
        Cs_ = c;
        *flags_.logstream << "DDC with nonzero S-const." << std::endl;
    }
    #endif
}

void DDE::initDDCConstraint(const FlowField& u) {
    if (!baseflow_)
        std::cerr << "DDE::initConstraint: Base flow has not been created." << std::endl;

    // Calculate Ubaseyy_ and related quantities
    UbulkBase_ = Ubase_.mean();
    ChebyCoeff Ubasey = diff(Ubase_);
    Ubaseyy_ = diff(Ubasey);

    WbulkBase_ = Wbase_.mean();
    ChebyCoeff Wbasey = diff(Wbase_);
    Wbaseyy_ = diff(Wbasey);

    ChebyCoeff Tbasey = diff(Tbase_);
    Tbaseyy_ = diff(Tbasey);

    ChebyCoeff Sbasey = diff(Sbase_);
    Sbaseyy_ = diff(Sbasey);

    // Determine actual Ubulk and dPdx from initial data Ubase + u.
    Real Ly = b_ - a_;
    ChebyCoeff u00(My_, a_, b_, Spectral);
    Real reucmplx = 0;
    for (int ny = 0; ny < My_; ++ny) {
        if (u.taskid() == 0)
            reucmplx = Re(u.cmplx(0, ny, 0, 0));
#ifdef HAVE_MPI
        MPI_Bcast(&reucmplx, 1, MPI_DOUBLE, u.task_coeff(0, 0), *u.comm_world());
#endif
        u00[ny] = reucmplx;
    }
    ChebyCoeff du00dy = diff(u00);
    UbulkAct_ = UbulkBase_ + u00.mean();
    dPdxAct_ = flags_.nu * (du00dy.eval_b() - du00dy.eval_a()) / Ly;
    ChebyCoeff w00(My_, a_, b_, Spectral);
    for (int ny = 0; ny < My_; ++ny) {
        if (u.taskid() == 0)
            reucmplx = Re(u.cmplx(0, ny, 0, 2));
#ifdef HAVE_MPI
        MPI_Bcast(&reucmplx, 1, MPI_DOUBLE, u.task_coeff(0, 0), *u.comm_world());
#endif
        w00[ny] = reucmplx;
    }
    ChebyCoeff dw00dy = diff(w00);
    WbulkAct_ = WbulkBase_ + w00.mean();
    dPdzAct_ = flags_.nu * (dw00dy.eval_b() - dw00dy.eval_a()) / Ly;

#ifdef HAVE_MPI
    MPI_Bcast(&dPdxAct_, 1, MPI_DOUBLE, u.task_coeff(0, 0), *u.comm_world());
    MPI_Bcast(&UbulkAct_, 1, MPI_DOUBLE, u.task_coeff(0, 0), *u.comm_world());
    MPI_Bcast(&dPdzAct_, 1, MPI_DOUBLE, u.task_coeff(0, 0), *u.comm_world());
    MPI_Bcast(&WbulkAct_, 1, MPI_DOUBLE, u.task_coeff(0, 0), *u.comm_world());
#endif

    if (Ubase_.length() != 0)
        dPdxAct_ += flags_.nu * (Ubasey.eval_b() - Ubasey.eval_a()) / Ly;
    if (Wbase_.length() != 0)
        dPdzAct_ += flags_.nu * (Wbasey.eval_b() - Wbasey.eval_a()) / Ly;

    if (flags_.constraint == BulkVelocity) {
        UbulkRef_ = flags_.Ubulk;
        WbulkRef_ = flags_.Wbulk;
    } else {
        dPdxAct_ = flags_.dPdx;
        dPdxRef_ = flags_.dPdx;
        dPdzAct_ = flags_.dPdz;
        dPdzRef_ = flags_.dPdz;
    }

    constraint_ = true;
}


ChebyCoeff laminarVelocityProfile(Real gammax, Real dPdx, Real Ubulk, Real Ua, Real Ub, Real a, Real b, int Ny,
                                  DDCFlags flags) {
    MeanConstraint constraint = flags.constraint;
    Real Vsuck = flags.Vsuck;

    Real Rey = flags.Rey;
    Real Pr = flags.Pr;
    Real Ra = flags.Ra;
    Real Rrho = flags.Rrho;
    Real Rsep = flags.Rsep;
    Real Ri = flags.Ri;
    Real sgammax = sin(gammax);
    
    Real h = b-a;
    Real minor2 = b*b-a*a;
    Real minor3 = b*b*b-a*a*a;
    Real plus2 = b*b+a*a;
    Real plus3 = b*b*b+a*a*a;
    Real c = 0.5*(b-a);
    Real d = 0.5*(b+a);
    
    Real Ta = flags.tlowerwall;
    Real Tb = flags.tupperwall;
    Real Sa = flags.slowerwall;
    Real Sb = flags.supperwall;

    Real pt1 = (Tb-Ta)/h;
    Real pt0 = (b*Ta-a*Tb)/h;

    Real ps1 = (Sb-Sa)/h;
    Real ps0 = (b*Sa-a*Sb)/h;
    #if defined(P6)
    Real p3 = -1.0*P2*sgammax/(6.0*P1) * (P3*pt1-P4*ps1);
    Real p2 = -1.0*P2*sgammax/(2.0*P1) * (P3*pt0-P4*ps0);
    Real p1 = (Ub-Ua)/h 
            + P2*sgammax/(2.0*P1*h)*(1.0/3.0*(P3*pt1-P4*ps1)*minor3 + (P3*pt0-P4*ps0)*minor2);
    Real p0 = 0.5*((Ub+Ua)-(Ub-Ua)*(b+a)/h)
            + P2*sgammax/(12.0*P1)*(P3*pt1-P4*ps1)*(plus3-minor3*(b+a)/h)
            + P2*sgammax/(4.0*P1)*(P3*pt0-P4*ps0)*(plus2-minor2*(b+a)/h);
    #elif defined(P5)
    Real p3 = -1.0*P2*sgammax/(6.0*P1) * (P3*pt1);
    Real p2 = -1.0*P2*sgammax/(2.0*P1) * (P3*pt0);
    Real p1 = (Ub-Ua)/h 
            + P2*sgammax/(2.0*P1*h)*(1.0/3.0*(P3*pt1)*minor3 + (P3*pt0)*minor2);
    Real p0 = 0.5*((Ub+Ua)-(Ub-Ua)*(b+a)/h)
            + P2*sgammax/(12.0*P1)*(P3*pt1)*(plus3-minor3*(b+a)/h)
            + P2*sgammax/(4.0*P1)*(P3*pt0)*(plus2-minor2*(b+a)/h);
    #else
    Real p3 = 0;
    Real p2 = 0;
    Real p1 = (Ub-Ua)/h;
    Real p0 = 0.5*((Ub+Ua)-(Ub-Ua)*(b+a)/h);
    #endif
    // printf("%f*y^3+%f*y^2+%f*y+%f\n",p3,p2,p1,p0);fflush(stdout);

    ChebyCoeff u(Ny, a, b, Spectral);

    if (constraint == BulkVelocity) {
        cferror("Using DDC with constraint BulkVelocity is not implemented yet");
    } else {
        if (Vsuck < 1e-14) {
            if (dPdx < 1e-14) {
                u[0] =    (1.5*c*c*d+d*d*d)*p3+(0.5*c*c+d*d)*p2+d*p1+p0;// See documentation
                u[1] = (0.75*c*c*c+3*c*d*d)*p3        +2*c*d*p2+c*p1;
                u[2] =            1.5*c*c*d*p3      +0.5*c*c*p2;
                u[3] =           0.25*c*c*c*p3;
            } else {
                cferror("Using DDC with nonzero dPdx is not implemented yet");
            }
        } else {
            cferror("Using DDC with SuctionVelocity is not implemented yet");
        }
    }
    return u;
}
ChebyCoeff linearTemperatureProfile(Real a, Real b, int Ny, DDCFlags flags) {
    MeanConstraint constraint = flags.constraint;

    ChebyCoeff T(Ny, a, b, Spectral);
    Real Vsuck = flags.Vsuck;
    Real dPdx = flags.dPdx;
    Real Ta = flags.tlowerwall;
    Real Tb = flags.tupperwall;

    Real h = b-a;
    Real d = 0.5*(b+a);

    if (constraint == BulkVelocity) {
        cferror("Using DDC with constraint BulkVelocity is not implemented yet");
    } else {
        if (Vsuck < 1e-14) {
            if (dPdx < 1e-14) {
                T[0] = (d-a)/h * (Tb - Ta) + Ta;  // the boundary conditions are given in units of Delta_T=Ta-Tb
                T[1] = (d-a)/h * (Tb - Ta);
            } else {
                cferror("Using DDC with nonzero dPdx is not implemented yet");
            }
        } else {
            cferror("Using DDC with SuctionVelocity is not implemented yet");
        }
    }
    return T;
}

ChebyCoeff linearSalinityProfile(Real a, Real b, int Ny, DDCFlags flags) {
    MeanConstraint constraint = flags.constraint;

    ChebyCoeff S(Ny, a, b, Spectral);
    Real Vsuck = flags.Vsuck;
    Real dPdx = flags.dPdx;
    Real Sa = flags.slowerwall;
    Real Sb = flags.supperwall;

    Real h = b-a;
    Real d = 0.5*(b+a);

    if (constraint == BulkVelocity) {
        cferror("Using DDC with constraint BulkVelocity is not implemented yet");
    } else {
        if (Vsuck < 1e-14) {
            if (dPdx < 1e-14) {
                S[0] = (d-a)/h * (Sb - Sa) + Sa;  // the boundary conditions are given in units of Delta_T=Ta-Tb
                S[1] = (d-a)/h * (Sb - Sa);
            } else {
                cferror("Using DDC with nonzero dPdx is not implemented yet");
            }
        } else {
            cferror("Using DDC with SuctionVelocity is not implemented yet");
        }
    }
    return S;
}
ChebyCoeff hydrostaticPressureGradientY(ChebyCoeff Tbase, ChebyCoeff Sbase, DDCFlags flags) {
    Real Rey = flags.Rey;
    Real Pr = flags.Pr;
    Real Ra = flags.Ra;
    Real Rrho = flags.Rrho;
    Real Rsep = flags.Rsep;
    Real Ri = flags.Ri;

    MeanConstraint constraint = flags.constraint;
    Real Vsuck = flags.Vsuck;
    Real dPdx = flags.dPdx;

    Real cgamma = cos(flags.gammax);

    ChebyCoeff Py(Tbase);

    if (constraint == BulkVelocity) {
        cferror("Using DDC with constraint BulkVelocity is not implemented yet");
    } else {
        if (Vsuck < 1e-14) {
            if (dPdx < 1e-14) {
                // dxP(y) = P2*(P3*T0(y)-P4*S0(y))*cos(gammaX)
                #if defined(P6)
                Py[0] = P2*(P3*Tbase[0]-P4*Sbase[0])*cgamma;  // See documentation about base solution
                Py[1] = P2*(P3*Tbase[1]-P4*Sbase[1])*cgamma;
                #elif defined(P5)
                Py[0] = P2*P3*Tbase[0]*cgamma;  // See documentation about base solution
                Py[1] = P2*P3*Tbase[1]*cgamma;
                #endif
            } else {
                cferror("Using DDC with nonzero dPdx is not implemented yet");
            }
        } else {
            cferror("Using DDC with SuctionVelocity is not implemented yet");
        }
    }

    return Py;
}

}  // namespace chflow