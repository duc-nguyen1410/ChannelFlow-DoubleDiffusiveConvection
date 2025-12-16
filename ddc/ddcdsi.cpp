/**
 * Dynamical System Interface for the DDC module
 *
 * Original author: Duc Nguyen
 */


#include "modules/ddc/ddcdsi.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "cfbasics/mathdefs.h"
#include "channelflow/diffops.h"
// #include "viscoelastic/veutils.h"

#include "modules/ddc/ddc.h"

using namespace std;

namespace chflow {

/*utility functions*/

std::vector<Real> ddcstats(const FlowField& u, const FlowField& temp, const FlowField& salt, const DDCFlags flags) {
    double l2n = L2Norm(u);
    if (std::isnan(l2n)) {
        cferror("L2Norm(u) is nan");
    }

    FlowField u_tot = totalVelocity(u, flags);
    FlowField temp_tot = totalTemperature(temp, flags);
    FlowField salt_tot = totalSalinity(salt, flags);

    std::vector<Real> stats;
    Real KE = 0.5*L2Norm2(u); stats.push_back(KE);
    Real PE = 0.5*flags.Ri*L2Norm2(temp); stats.push_back(PE);
    stats.push_back(KE+PE);
    stats.push_back(dissipation(u_tot,flags));

    stats.push_back(wallshearUpper(u_tot)); // for Langham2019JFM
    stats.push_back(wallshear(u_tot)); // 
    stats.push_back(L2Norm(u));
    stats.push_back(L2Norm(u_tot));
    stats.push_back(L2Norm3d(u));
    stats.push_back(Ecf(u));
    stats.push_back(getUbulk(u));
    stats.push_back(getWbulk(u));

    stats.push_back(L2Norm(temp));
    stats.push_back(L2Norm(temp_tot));
    stats.push_back(heatcontent(temp_tot, flags));// averaged temp
    stats.push_back(heatflux(temp_tot, flags));// J_t
    stats.push_back(Nusselt_t_plane(u_tot, temp_tot, flags));// Nu_t

    stats.push_back(L2Norm(salt));
    stats.push_back(L2Norm(salt_tot));
    stats.push_back(saltcontent(salt_tot, flags));// averaged salt
    #ifdef P6
    stats.push_back(massflux(salt_tot, flags));// J_c
    stats.push_back(Nusselt_c_plane(u_tot, salt_tot, flags));// Nu_c
    #endif
    // stats.push_back(abs(temp_tot.dudy_a()));
    // stats.push_back(abs(salt_tot.dudy_a()));

    stats.push_back(buoyPowerInput(u_tot, temp_tot, salt_tot, flags));

    // stats.push_back(UdPcontent(u_tot, flags));
    
    return stats;
}

string ddcfieldstatsheader(const DDCFlags flags) {
    stringstream header;
    header  << setw(14) << "KinEnergy"
            << setw(14) << "PotEnergy"
            << setw(14) << "TotEnergy"
            << setw(14) << "Dissipation"

            << setw(14) << "wshearUpper" 
            << setw(14) << "wallshear" 
            << setw(14) << "L2(u')" 
            << setw(14) << "L2(u)" 
            << setw(14) << "e3d" 
            << setw(14) << "ecf" 
            << setw(14) << "ubulk" 
            << setw(14) << "wbulk" 

            << setw(14) << "L2(T')" 
            << setw(14) << "L2(T)" 
            << setw(14) << "<T>a" 
            << setw(14) << "heatflux" 
            << setw(14) << "Nu_t" 

            << setw(14) << "L2(S')" 
            << setw(14) << "L2(S)" 
            << setw(14) << "<S>a"
            #ifdef P6
            << setw(14) << "massflux" 
            << setw(14) << "Nu_c" 
            #endif
            // << setw(14) << "Nu" 
            // << setw(14) << "Sh"

            << setw(14) << "buoyPowIn"

            // << setw(14) << "UdP"
            ;
    return header.str();
}

string ddcfieldstatsheader_t(const string tname, const DDCFlags flags) {
    stringstream header;
    header << setw(6) << "#(" << tname << ")" << ddcfieldstatsheader(flags);
    return header.str();
}

string ddcfieldstats(const FlowField& u, const FlowField& temp, const FlowField& salt, const DDCFlags flags) {
    std::vector<Real> stats = ddcstats(u, temp, salt, flags);
    // Return string
    stringstream s;
    for (uint i = 0; i < stats.size(); i++) {
        s << setw(14) << stats[i];
    }
    return s.str();
}

string ddcfieldstats_t(const FlowField& u, const FlowField& temp, const FlowField& salt, const Real t, const DDCFlags flags) {
    std::vector<Real> stats = ddcstats(u, temp, salt, flags);
    // Return string
    stringstream s;
    s << setw(8) << t;
    for (uint i = 0; i < stats.size(); i++) {
        s << setw(14) << stats[i];
    }
    return s.str();
}

FlowField totalVelocity(const FlowField& velo, const DDCFlags flags) {
    // copy
    FlowField u(velo);
    FlowField tmp(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());

    // get base flow
    DDC ddc({u, tmp, tmp, tmp}, flags);
    ChebyCoeff Ubase = ddc.Ubase();
    ChebyCoeff Wbase = ddc.Wbase();

    // add base flow (should be identical to code in function temperatureNL in DDE
    for (int ny = 0; ny < u.Ny(); ++ny) {
        if (u.taskid() == u.task_coeff(0, 0)) {
            u.cmplx(0, ny, 0, 0) += Complex(Ubase(ny), 0.0);
            u.cmplx(0, ny, 0, 2) += Complex(Wbase(ny), 0.0);
        }
    }
    if (u.taskid() == u.task_coeff(0, 0)) {
        u.cmplx(0, 0, 0, 1) -= Complex(flags.Vsuck, 0.);
    }

    return u;
}

FlowField totalTemperature(const FlowField& temp, const DDCFlags flags) {
    // copy
    FlowField T(temp);
    FlowField tmp(T.Nx(), T.Ny(), T.Nz(), 3, T.Lx(), T.Lz(), T.a(), T.b(), T.cfmpi());

    // get base flow
    DDC ddc({tmp, T, T, T}, flags);
    ChebyCoeff Tbase = ddc.Tbase();

    // add base flow (should be identical to code in function temperatureNL in DDE
    for (int ny = 0; ny < T.Ny(); ++ny) {
        if (T.taskid() == T.task_coeff(0, 0))
            T.cmplx(0, ny, 0, 0) += Complex(Tbase(ny), 0.0);
    }

    return T;
}

FlowField totalSalinity(const FlowField& salt, const DDCFlags flags) {
    // copy
    FlowField S(salt);
    FlowField tmp(S.Nx(), S.Ny(), S.Nz(), 3, S.Lx(), S.Lz(), S.a(), S.b(), S.cfmpi());

    // get base flow
    DDC ddc({tmp, S, S, S}, flags);
    ChebyCoeff Sbase = ddc.Sbase();

    // add base flow (should be identical to code in function temperatureNL in DDE
    for (int ny = 0; ny < S.Ny(); ++ny) {
        if (S.taskid() == S.task_coeff(0, 0))
            S.cmplx(0, ny, 0, 0) += Complex(Sbase(ny), 0.0);
    }

    return S;
}


Real heatcontent(const FlowField& ttot, const DDCFlags flags) {
    assert(ttot.ystate() == Spectral);
    int N = 100;
    Real dy = (flags.ystats - ttot.a()) / (N - 1);
    Real avt = 0;  // average temperature
    if (ttot.taskid() == ttot.task_coeff(0, 0)) {
        ChebyCoeff tprof = Re(ttot.profile(0, 0, 0));
        for (int i = 0; i < N; ++i) {
            Real y = ttot.a() + i * dy;
            avt += tprof.eval(y);
        }
        avt *= 1.0 / N;
    }
#ifdef HAVE_MPI
    MPI_Bcast(&avt, 1, MPI_DOUBLE, ttot.task_coeff(0, 0), ttot.cfmpi()->comm_world);
#endif
    return avt;
}
Real saltcontent(const FlowField& stot, const DDCFlags flags) {
    assert(stot.ystate() == Spectral);
    int N = 100;
    Real dy = (flags.ystats - stot.a()) / (N - 1);
    Real avt = 0;  // average salinity
    if (stot.taskid() == stot.task_coeff(0, 0)) {
        ChebyCoeff sprof = Re(stot.profile(0, 0, 0));
        for (int i = 0; i < N; ++i) {
            Real y = stot.a() + i * dy;
            avt += sprof.eval(y);
        }
        avt *= 1.0 / N;
    }
#ifdef HAVE_MPI
    MPI_Bcast(&avt, 1, MPI_DOUBLE, stot.task_coeff(0, 0), stot.cfmpi()->comm_world);
#endif
    return avt;
}

// Real UdPcontent(const FlowField& utot, const DDCFlags flags) {

//     FlowField u(utot);
//     Vector y = u.ygridpts();
//     // FlowField dP(Pbasey_);
//     // FlowField dP(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Physical, Physical);
//     // dP+=flags.Pbasey_;
//     FlowField yinput(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Physical, Physical);
//     lint Nz = u.Nz();
//     lint nxlocmin = u.nxlocmin();
//     lint nxlocmax = u.nxlocmin() + u.Nxloc();
//     lint nylocmin = u.nylocmin();
//     lint nylocmax = u.nylocmax();

//     u.makePhysical();
//     // dP.makePhysical();
// #ifdef HAVE_MPI
//     for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
//         for (lint nz = 0; nz < Nz; ++nz)
//             for (lint ny = nylocmin; ny < nylocmax; ++ny) {
//                 yinput(nx, ny, nz, 0) = u(nx, ny, nz, 1) * y(ny);
//             }
// #else
//     for (lint ny = nylocmin; ny < nylocmax; ++ny)
//         for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
//             for (lint nz = 0; nz < Nz; ++nz) {
//                 yinput(nx, ny, nz, 0) = u(nx, ny, nz, 1) * y(ny);
//             }
// #endif
//     yinput.makeSpectral();

//     ChebyCoeff yprof(u.My(), u.a(), u.b(), Spectral);
//     Real ytmp = 0;
//     for (int ny = 0; ny < u.My(); ++ny) {
//         if (u.taskid() == 0) {
//             ytmp = Re(yinput.cmplx(0, ny, 0, 0)); 
//         }
// #ifdef HAVE_MPI
//         MPI_Bcast(&ytmp, 1, MPI_DOUBLE, u.task_coeff(0, 0), *u.comm_world());
// #endif
//         yprof[ny] = ytmp;
//     }

//     Real output = yprof.mean();
// #ifdef HAVE_MPI
//     MPI_Bcast(&output, 1, MPI_DOUBLE, u.task_coeff(0, 0), u.cfmpi()->comm_world);
// #endif
//     return output;

//     // assert(utot.ystate() == Spectral);

//     // ChebyCoeff uprof = Re(utot.profile(0, 0, 0));
//     // ChebyCoeff UdP = uprof[1]*Py;

//     // int N = 100;
//     // Real dy = (flags.ystats - stot.a()) / (N - 1);
//     // Real avt = 0;  // average salinity
//     // if (stot.taskid() == stot.task_coeff(0, 0)) {
//     //     ChebyCoeff sprof = Re(stot.profile(0, 0, 0));
//     //     for (int i = 0; i < N; ++i) {
//     //         Real y = stot.a() + i * dy;
//     //         avt += sprof.eval(y);
//     //     }
//     //     avt *= 1.0 / N;
//     // }
// // #ifdef HAVE_MPI
// //     MPI_Bcast(&avt, 1, MPI_DOUBLE, stot.task_coeff(0, 0), stot.cfmpi()->comm_world);
// // #endif
// //     return avt;
// }

Real dissipation(const FlowField& utot, const DDCFlags flags, bool normalize, bool relative) {
    // get parameters
    Real Rey = flags.Rey;
    Real Pr = flags.Pr;
    Real Ra = flags.Ra;
    Real Rrho = flags.Rrho;
    Real Rsep = flags.Rsep;
    Real Ri = flags.Ri;

    Real nu = P1;
    
    Real diss = nu * dissipation(utot, normalize);
    // printf("dissipation=%f, nu=%f, diss=%f\n",dissipation(utot, normalize),flags.nu,diss);fflush(stdout);
    // analytic laminar dissipation (for standard base flow)
    Real sing = sin(flags.gammax);
    Real grav = 1.0;
    Real laminarDiss = grav * sing * sing / (720 * nu);  // normalized by Volume
    // difference between full input and laminar input
    if (relative && abs(laminarDiss) > 1e-12) {
        diss *= 1.0 / laminarDiss;
        diss -= 1;
    }

    return diss;
}

Real heatflux(const FlowField& temp, const DDCFlags flags, bool normalize, bool relative) {
    // with reference to wallshear, but only lower wall
    assert(temp.ystate() == Spectral);

    // get parameters
    Real Rey = flags.Rey;
    Real Pr = flags.Pr;
    Real Ra = flags.Ra;
    Real Rrho = flags.Rrho;
    Real Rsep = flags.Rsep;
    Real Ri = flags.Ri;

    Real kappa = P5;
    Real I = 0;
    if (temp.taskid() == temp.task_coeff(0, 0)) {
        ChebyCoeff tprof = Re(temp.profile(0, 0, 0));
        ChebyCoeff dTdy = diff(tprof);
        I = kappa * abs(dTdy.eval_a());
    }
#ifdef HAVE_MPI
    MPI_Bcast(&I, 1, MPI_DOUBLE, temp.task_coeff(0, 0), temp.cfmpi()->comm_world);
#endif
    if (!normalize)
        I *= 2 * temp.Lx() * temp.Lz();
    if (relative) {
        Real deltaT = flags.tlowerwall - flags.tupperwall;
        Real H = temp.b() - temp.a();
        I *= 1.0 / kappa / deltaT * H;
        // I -= 1.0;
    }
    return I;
}
#ifdef P6
Real massflux(const FlowField& temp, const DDCFlags flags, bool normalize, bool relative) {
    // with reference to wallshear, but only lower wall
    assert(temp.ystate() == Spectral);

    // get parameters
    Real Rey = flags.Rey;
    Real Pr = flags.Pr;
    Real Ra = flags.Ra;
    Real Rrho = flags.Rrho;
    Real Rsep = flags.Rsep;
    Real Ri = flags.Ri;
    Real Le = flags.Le;

    Real kappa = P6;
    Real I = 0;
    if (temp.taskid() == temp.task_coeff(0, 0)) {
        ChebyCoeff tprof = Re(temp.profile(0, 0, 0));
        ChebyCoeff dTdy = diff(tprof);
        I = kappa * abs(dTdy.eval_a());
    }
#ifdef HAVE_MPI
    MPI_Bcast(&I, 1, MPI_DOUBLE, temp.task_coeff(0, 0), temp.cfmpi()->comm_world);
#endif
    if (!normalize)
        I *= 2 * temp.Lx() * temp.Lz();
    if (relative) {
        Real deltaC = flags.slowerwall - flags.supperwall;
        Real H = temp.b() - temp.a();
        I *= 1.0 / kappa / deltaC * H;
        // I -= 1.0;
    }
    return I;
}
#endif
Real Nusselt_t_plane(const FlowField& utot, const FlowField& ttot, const DDCFlags flags, bool relative) {
    assert(utot.ystate() == Spectral && ttot.ystate() == Spectral);

    // get parameters
    Real Rey = flags.Rey;
    Real Pr = flags.Pr;
    Real Ra = flags.Ra;
    Real Rrho = flags.Rrho;
    Real Rsep = flags.Rsep;
    Real Ri = flags.Ri;

    // calculate product for advective heat transport
    FlowField u(utot);
    FlowField T(ttot);
    FlowField vt(T.Nx(), T.Ny(), T.Nz(), T.Nd(), T.Lx(), T.Lz(), T.a(), T.b(), T.cfmpi(), Physical, Physical);
    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    // loop to form product
    u.makePhysical();
    T.makePhysical();
#ifdef HAVE_MPI
    for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
        for (lint nz = 0; nz < Nz; ++nz)
            for (lint ny = nylocmin; ny < nylocmax; ++ny) {
                vt(nx, ny, nz, 0) = u(nx, ny, nz, 1) * T(nx, ny, nz, 0);
            }
#else
    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (lint nz = 0; nz < Nz; ++nz) {
                vt(nx, ny, nz, 0) = u(nx, ny, nz, 1) * T(nx, ny, nz, 0);
            }
#endif
    vt.makeSpectral();

    // calculate Nusselt number
    Real Nu = 0;
    Real y = flags.ystats;
    Real kappa = P5;
    if (ttot.taskid() == ttot.task_coeff(0, 0)) {
        ChebyCoeff tprof = Re(ttot.profile(0, 0, 0));
        ChebyCoeff dTdy = diff(tprof);
        ChebyCoeff vtprof = Re(vt.profile(0, 0, 0));
        Nu = vtprof.eval(y) -
             kappa * dTdy.eval(y);  // Formula 10 in Chilla&Schumacher 2012 (together with normalization below)
    }

#ifdef HAVE_MPI
    MPI_Bcast(&Nu, 1, MPI_DOUBLE, ttot.task_coeff(0, 0), ttot.cfmpi()->comm_world);
#endif

    if (relative) {
        Real deltaT = flags.tlowerwall - flags.tupperwall;
        Real H = utot.b() - utot.a();
        Nu *= 1.0 / kappa / deltaT * H;
        // Nu -= 1.0;
    }
    return Nu;
}

#ifdef P6
Real Nusselt_c_plane(const FlowField& utot, const FlowField& ttot, const DDCFlags flags, bool relative) {
    assert(utot.ystate() == Spectral && ttot.ystate() == Spectral);

    // get parameters
    Real Rey = flags.Rey;
    Real Pr = flags.Pr;
    Real Ra = flags.Ra;
    Real Rrho = flags.Rrho;
    Real Rsep = flags.Rsep;
    Real Ri = flags.Ri;
    Real Le = flags.Le;

    // calculate product for advective heat transport
    FlowField u(utot);
    FlowField T(ttot);
    FlowField vt(T.Nx(), T.Ny(), T.Nz(), T.Nd(), T.Lx(), T.Lz(), T.a(), T.b(), T.cfmpi(), Physical, Physical);
    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    // loop to form product
    u.makePhysical();
    T.makePhysical();
#ifdef HAVE_MPI
    for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
        for (lint nz = 0; nz < Nz; ++nz)
            for (lint ny = nylocmin; ny < nylocmax; ++ny) {
                vt(nx, ny, nz, 0) = u(nx, ny, nz, 1) * T(nx, ny, nz, 0);
            }
#else
    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (lint nz = 0; nz < Nz; ++nz) {
                vt(nx, ny, nz, 0) = u(nx, ny, nz, 1) * T(nx, ny, nz, 0);
            }
#endif
    vt.makeSpectral();

    // calculate Nusselt number
    Real Nu = 0;
    Real y = flags.ystats;
    Real kappa = P6;
    if (ttot.taskid() == ttot.task_coeff(0, 0)) {
        ChebyCoeff tprof = Re(ttot.profile(0, 0, 0));
        ChebyCoeff dTdy = diff(tprof);
        ChebyCoeff vtprof = Re(vt.profile(0, 0, 0));
        Nu = vtprof.eval(y) -
             kappa * dTdy.eval(y);  // Formula 10 in Chilla&Schumacher 2012 (together with normalization below)
    }

#ifdef HAVE_MPI
    MPI_Bcast(&Nu, 1, MPI_DOUBLE, ttot.task_coeff(0, 0), ttot.cfmpi()->comm_world);
#endif

    if (relative) {
        Real deltaC = flags.slowerwall - flags.supperwall;
        Real H = utot.b() - utot.a();
        Nu *= 1.0 / kappa / deltaC * H;
        // I -= 1.0;
    }
    return Nu;
}
#endif

Real buoyPowerInput(const FlowField& utot, const FlowField& ttot, const FlowField& stot, const DDCFlags flags, bool relative) {
    // calculate the bouyancy force along the velocity field to get
    // the power input to the kinetic energy equation

    // get parameters
    Real Rey = flags.Rey;
    Real Pr = flags.Pr;
    Real Ra = flags.Ra;
    Real Rrho = flags.Rrho;
    Real Rsep = flags.Rsep;
    Real Ri = flags.Ri;
    
    Real nu = P1;

    Real sing = sin(flags.gammax);
    Real cosg = cos(flags.gammax);
    Real grav = 1.0; 
    Real laminarInput = grav * sing * sing / (720 * nu);  // normalized by Volume

    // prepare loop over field
    FlowField u(utot);
    FlowField T(ttot);
    FlowField S(stot);
    FlowField xinput(T.Nx(), T.Ny(), T.Nz(), T.Nd(), T.Lx(), T.Lz(), T.a(), T.b(), T.cfmpi(), Physical, Physical);
    FlowField yinput(T.Nx(), T.Ny(), T.Nz(), T.Nd(), T.Lx(), T.Lz(), T.a(), T.b(), T.cfmpi(), Physical, Physical);
    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    // sum up buoyancy term U*F = p2(p3*T-p4*S)*U
    u.makePhysical();
    T.makePhysical();
    S.makePhysical();
#ifdef HAVE_MPI
    for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
        for (lint nz = 0; nz < Nz; ++nz)
            for (lint ny = nylocmin; ny < nylocmax; ++ny) {
                #ifdef P6
                xinput(nx, ny, nz, 0) = u(nx, ny, nz, 0) * P2*(P3*T(nx, ny, nz, 0) - P4*S(nx, ny, nz, 0));
                yinput(nx, ny, nz, 0) = u(nx, ny, nz, 1) * P2*(P3*T(nx, ny, nz, 0) - P4*S(nx, ny, nz, 0));
                #elif defined(P5)
                xinput(nx, ny, nz, 0) = u(nx, ny, nz, 0) * P2*P3*T(nx, ny, nz, 0);
                yinput(nx, ny, nz, 0) = u(nx, ny, nz, 1) * P2*P3*T(nx, ny, nz, 0);
                #endif
            }
#else
    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (lint nz = 0; nz < Nz; ++nz) {
                #ifdef P6
                xinput(nx, ny, nz, 0) = u(nx, ny, nz, 0) * P2*(P3*T(nx, ny, nz, 0) - P4*S(nx, ny, nz, 0));
                yinput(nx, ny, nz, 0) = u(nx, ny, nz, 1) * P2*(P3*T(nx, ny, nz, 0) - P4*S(nx, ny, nz, 0));
                #elif defined(P5)
                xinput(nx, ny, nz, 0) = u(nx, ny, nz, 0) * P2*P3*T(nx, ny, nz, 0);
                yinput(nx, ny, nz, 0) = u(nx, ny, nz, 1) * P2*P3*T(nx, ny, nz, 0);
                #endif
            }
#endif
    xinput.makeSpectral();
    yinput.makeSpectral();

    // calculate the input mean with cheby profile (code is taken from OBE::initConstraint)
    ChebyCoeff xprof(T.My(), T.a(), T.b(), Spectral);
    ChebyCoeff yprof(T.My(), T.a(), T.b(), Spectral);
    Real xtmp = 0;
    Real ytmp = 0;
    for (int ny = 0; ny < T.My(); ++ny) {
        if (utot.taskid() == 0) {
            xtmp = sing * Re(xinput.cmplx(0, ny, 0, 0)); // U*F * sin(gammax)*ex
            ytmp = cosg * Re(yinput.cmplx(0, ny, 0, 0)); // U*F * cos(gammax)*ey
        }
#ifdef HAVE_MPI
        MPI_Bcast(&xtmp, 1, MPI_DOUBLE, utot.task_coeff(0, 0), *utot.comm_world());
        MPI_Bcast(&ytmp, 1, MPI_DOUBLE, utot.task_coeff(0, 0), *utot.comm_world());
#endif
        xprof[ny] = xtmp;
        yprof[ny] = ytmp;
    }

    Real buoyancyInput = xprof.mean() + yprof.mean();
    // printf("xprof.mean()=%f, yprof.mean()=%f, buoyancyInput=%f\n",xprof.mean(),yprof.mean(),buoyancyInput);fflush(stdout);
    // difference between full input and laminar input
    if (relative && abs(laminarInput) > 1e-12) {
        buoyancyInput *= 1.0 / laminarInput;
        buoyancyInput -= 1;
    }

    return buoyancyInput;
}


/* Begin of ddcDSI class*/

ddcDSI::ddcDSI() {}

ddcDSI::ddcDSI(DDCFlags& ddcflags, FieldSymmetry sigma, PoincareCondition* h, TimeStep dt, bool Tsearch, bool xrelative,
               bool zrelative, bool Tnormalize, Real Unormalize, const FlowField& u, const FlowField& temp, const FlowField& salt, ostream* os)
    : cfDSI(ddcflags, sigma, h, dt, Tsearch, xrelative, zrelative, Tnormalize, Unormalize, u, os),
      ddcflags_(ddcflags) {}

Eigen::VectorXd ddcDSI::eval(const Eigen::VectorXd& x) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);

    Real T;
    extractVectorDDC(x, u, temp, salt, sigma_, T);

    FlowField Gu(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField Gtemp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField Gsalt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    G(u, temp, salt, T, h_, sigma_, Gu, Gtemp, Gsalt, ddcflags_, dt_, Tnormalize_, Unormalize_, fcount_, CFL_, *os_);
    Eigen::VectorXd Gx(Eigen::VectorXd::Zero(x.rows()));
    //   Galpha *= 1./vednsflags_.b_para;
    field2vector(Gu, Gtemp, Gsalt, Gx);  // This does not change the size of Gx and automatically leaves the last entries zero

    return Gx;
}

Eigen::VectorXd ddcDSI::eval(const Eigen::VectorXd& x0, const Eigen::VectorXd& x1, bool symopt) {
    FlowField u0(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField u1(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp0(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp1(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt0(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt1(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    Real T0, T1;
    FieldSymmetry sigma0, sigma1;
    extractVectorDDC(x0, u0, temp0, salt0, sigma0, T0);
    extractVectorDDC(x1, u1, temp1, salt1, sigma1, T1);

    FlowField Gu(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField Gtemp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField Gsalt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);

    f(u0, temp0, salt0, T0, h_, Gu, Gtemp, Gsalt, ddcflags_, dt_, fcount_, CFL_, *os_);
    if (symopt) {
        Gu *= sigma0;
        if (sigma0.sy() == -1) {
            // wall-normal mirroring in velocity requires sign change in temperature
            FieldSymmetry inv(-1);
            sigma0 *= inv;
        }
        Gtemp *= sigma0;
        Gsalt *= sigma0;
    }
    Gu -= u1;
    Gtemp -= temp1;
    Gsalt -= salt1;

    // normalize
    if (Tnormalize_) {
        Gu *= 1.0 / T0;
        Gtemp *= 1.0 / T0;
        Gsalt *= 1.0 / T0;
    }
    if (Unormalize_ != 0.0) {
        Real funorm = L2Norm3d(Gu);
        Gu *= 1. / sqrt(abs(funorm * (Unormalize_ - funorm)));
        // u should stay off zero, so normalize with u for now - temp should also stay away from zero
        Gtemp *= 1. / sqrt(abs(funorm * (Unormalize_ - funorm)));
        Gsalt *= 1. / sqrt(abs(funorm * (Unormalize_ - funorm)));
    }

    Eigen::VectorXd Gx(Eigen::VectorXd::Zero(x0.rows()));
    field2vector(Gu, Gtemp, Gsalt, Gx);  // This does not change the size of Gx and automatically leaves the last entries zero

    return Gx;
}

void ddcDSI::save(const Eigen::VectorXd& x, const string filebase, const string outdir, const bool fieldsonly) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FieldSymmetry sigma;
    Real T;
    extractVectorDDC(x, u, temp, salt, sigma, T);

    u.save(outdir + "u" + filebase);
    temp.save(outdir + "t" + filebase);
    salt.save(outdir + "s" + filebase);

    if (!fieldsonly) {
        string fs = ddcfieldstats(u, temp, salt, ddcflags_);
        if (u.taskid() == 0) {
            if (xrelative_ || zrelative_ || !sigma.isIdentity())
                sigma.save(outdir + "sigma" + filebase);
            if (Tsearch_)
                chflow::save(T, outdir + "T" + filebase);
            // sigma.save (outdir+"sigmaconverge.asc", ios::app);
            ofstream fout((outdir + "fieldconverge.asc").c_str(), ios::app);
            long pos = fout.tellp();
            if (pos == 0)
                fout << ddcfieldstatsheader() << endl;
            fout << fs << endl;
            fout.close();
            ddcflags_.save(outdir);
        }
    }
}

string ddcDSI::stats(const Eigen::VectorXd& x) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FieldSymmetry sigma;
    Real T;
    extractVectorDDC(x, u, temp, salt, sigma, T);
    return ddcfieldstats_t(u, temp, salt, mu_, ddcflags_);
}

pair<string, string> ddcDSI::stats_minmax(const Eigen::VectorXd& x) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField Gu(u);
    FlowField Gtemp(temp);
    FlowField Gsalt(salt);
    FieldSymmetry sigma;
    Real T;
    extractVectorDDC(x, u, temp, salt, sigma, T);

    std::vector<Real> stats = ddcstats(u, temp, salt, ddcflags_);
    std::vector<Real> minstats(stats);
    std::vector<Real> maxstats(stats);

    // quick hack to avoid new interface or creating simple f() for DDC
    TimeStep dt = TimeStep(ddcflags_.dt, 0, 1, 1, 0, 0, false);
    int fcount = 0;
    PoincareCondition* h = 0;
    Real CFL = 0.0;
    std::ostream muted_os(0);
    Real timep = T / 100.0;

    *os_ << "Using flag -orbOut: Calculate minmax-statistics of periodic orbit." << endl;
    for (int t = 0; t < 100; t++) {
        f(u, temp, salt, timep, h, Gu, Gtemp, Gsalt, ddcflags_, dt, fcount, CFL, muted_os);
        stats = ddcstats(Gu, Gtemp, Gsalt, ddcflags_);
        for (uint i = 0; i < stats.size(); i++) {
            minstats[i] = (minstats[i] < stats[i]) ? minstats[i] : stats[i];
            maxstats[i] = (maxstats[i] > stats[i]) ? maxstats[i] : stats[i];
        }
        u = Gu;
        temp = Gtemp;
        salt = Gsalt;
    }
    // Return string
    stringstream smin;
    stringstream smax;
    smin << setw(8) << mu_;
    smax << setw(8) << mu_;
    for (uint i = 0; i < stats.size(); i++) {
        smin << setw(14) << minstats[i];
        smax << setw(14) << maxstats[i];
    }

    pair<string, string> minmax;
    minmax = make_pair(smin.str(), smax.str());
    return minmax;
}

string ddcDSI::statsHeader() { return ddcfieldstatsheader_t(ddc_cPar2s(ddc_cPar_), ddcflags_); }

/// after finding new solution fix phases
void ddcDSI::phaseShift(Eigen::VectorXd& x) {
    if (xphasehack_ || zphasehack_) {
        FlowField unew(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
        FlowField tnew(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
        FlowField snew(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
        FieldSymmetry sigma;
        Real T;
        extractVectorDDC(x, unew, tnew, snew, sigma, T);
        // vector2field (x,unew);
        const int phasehackcoord = 0;  // Those values were fixed in continuesoln anyway
        const parity phasehackparity = Odd;
        const Real phasehackguess = 0.0;

        if (zphasehack_) {
            FieldSymmetry tau = zfixphasehack(unew, phasehackguess, phasehackcoord, phasehackparity);
            cout << "fixing z phase of potential solution with phase shift tau == " << tau << endl;
            unew *= tau;
            tnew *= tau;
            snew *= tau;
        }
        if (xphasehack_) {
            FieldSymmetry tau = xfixphasehack(unew, phasehackguess, phasehackcoord, phasehackparity);
            cout << "fixing x phase of potential solution with phase shift tau == " << tau << endl;
            unew *= tau;
            tnew *= tau;
            snew *= tau;
        }
        if (uUbasehack_) {
            cout << "fixing u+Ubase decomposition so that <du/dy> = 0 at walls (i.e. Ubase balances mean pressure "
                    "gradient))"
                 << endl;
            Real ubulk = Re(unew.profile(0, 0, 0)).mean();
            if (abs(ubulk) < 1e-15)
                ubulk = 0.0;

            ChebyCoeff Ubase = laminarProfile(ddcflags_.nu, ddcflags_.constraint, ddcflags_.dPdx,
                                              ddcflags_.Ubulk - ubulk, ddcflags_.Vsuck, unew.a(), unew.b(),
                                              ddcflags_.ulowerwall, ddcflags_.uupperwall, unew.Ny());

            fixuUbasehack(unew, Ubase);
        }
        makeVectorDDC(unew, tnew, snew, sigma, T, x);
    }
}

void ddcDSI::phaseShift(Eigen::MatrixXd& y) {
    if (xphasehack_ || zphasehack_) {
        FlowField unew(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
        FlowField tnew(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
        FlowField snew(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
        Eigen::VectorXd yvec;
        FieldSymmetry sigma;
        Real T;

        const int phasehackcoord = 0;  // Those values were fixed in continuesoln anyway
        const parity phasehackparity = Odd;
        const Real phasehackguess = 0.0;

        FieldSymmetry taux(0.0, 0.0);
        FieldSymmetry tauz(0.0, 0.0);

        extractVectorDDC(y.col(0), unew, tnew, snew, sigma, T);

        if (xphasehack_) {
            taux = xfixphasehack(unew, phasehackguess, phasehackcoord, phasehackparity);
            cout << "fixing x phase of potential solution with phase shift tau == " << taux << endl;
        }
        if (zphasehack_) {
            tauz = zfixphasehack(unew, phasehackguess, phasehackcoord, phasehackparity);
            cout << "fixing z phase of potential solution with phase shift tau == " << tauz << endl;
        }

        for (int i = 0; i < y.cols(); i++) {
            extractVectorDDC(y.col(i), unew, tnew, snew, sigma, T);
            unew *= taux;
            tnew *= taux;
            snew *= taux;
            unew *= tauz;
            tnew *= tauz;
            snew *= tauz;
            makeVectorDDC(unew, tnew, snew, sigma, T, yvec);
            y.col(i) = yvec;
        }
    }
}

Real ddcDSI::extractT(const Eigen::VectorXd& x) {  // inefficient hack
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FieldSymmetry sigma;
    Real T;
    extractVectorDDC(x, u, temp, salt, sigma, T);
    return T;
}

Real ddcDSI::extractXshift(const Eigen::VectorXd& x) {  // inefficient hack
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FieldSymmetry sigma;
    Real T;
    extractVectorDDC(x, u, temp, salt, sigma, T);
    return sigma.ax();
}

Real ddcDSI::extractZshift(const Eigen::VectorXd& x) {  // inefficient hack
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FieldSymmetry sigma;
    Real T;
    extractVectorDDC(x, u, temp, salt, sigma, T);
    return sigma.az();
}

void ddcDSI::makeVectorDDC(const FlowField& u, const FlowField& temp, const FlowField& salt, const FieldSymmetry& sigma, const Real T,
                           Eigen::VectorXd& x) {
    if (u.Nd() != 3)
        cferror("ddcDSI::makeVector(): u.Nd() = " + i2s(u.Nd()) + " != 3");
    if (temp.Nd() != 1)
        cferror("ddcDSI::makeVector(): temp.Nd() = " + i2s(temp.Nd()) + " != 1");
    if (salt.Nd() != 1)
        cferror("ddcDSI::makeVector(): salt.Nd() = " + i2s(salt.Nd()) + " != 1");
    int taskid = u.taskid();

    int uunk = field2vector_size(u, temp, salt);                   // # of variables for u and alpha unknonwn
    const int Tunk = (Tsearch_ && taskid == 0) ? uunk : -1;  // index for T unknown
    const int xunk = (xrelative_ && taskid == 0) ? uunk + Tsearch_ : -1;
    const int zunk = (zrelative_ && taskid == 0) ? uunk + Tsearch_ + xrelative_ : -1;
    int Nunk = (taskid == 0) ? uunk + Tsearch_ + xrelative_ + zrelative_ : uunk;
    if (x.rows() < Nunk)
        x.resize(Nunk);
    field2vector(u, temp, salt, x);
    if (taskid == 0) {
        if (Tsearch_)
            x(Tunk) = T;
        if (xrelative_)
            x(xunk) = sigma.ax();
        if (zrelative_)
            x(zunk) = sigma.az();
    }
}

void ddcDSI::extractVectorDDC(const Eigen::VectorXd& x, FlowField& u, FlowField& temp, FlowField& salt, FieldSymmetry& sigma, Real& T) {
    int uunk = field2vector_size(u, temp, salt);  // number of components in x that corresond to u and alpha
    vector2field(x, u, temp, salt);
    const int Tunk = uunk + Tsearch_ - 1;
    const int xunk = uunk + Tsearch_ + xrelative_ - 1;
    const int zunk = uunk + Tsearch_ + xrelative_ + zrelative_ - 1;
    Real ax = 0;
    Real az = 0;
    if (u.taskid() == 0) {
        T = Tsearch_ ? x(Tunk) : Tinit_;
        ax = xrelative_ ? x(xunk) : axinit_;
        az = zrelative_ ? x(zunk) : azinit_;
    }
#ifdef HAVE_MPI
    MPI_Bcast(&T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&az, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    sigma = FieldSymmetry(sigma_.sx(), sigma_.sy(), sigma_.sz(), ax, az, sigma_.s());
}

Eigen::VectorXd ddcDSI::xdiff(const Eigen::VectorXd& a) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    vector2field(a, u, temp, salt);
    Eigen::VectorXd dadx(a.size());
    dadx.setZero();
    u = chflow::xdiff(u);
    temp = chflow::xdiff(temp);
    salt = chflow::xdiff(salt);
    field2vector(u, temp, salt, dadx);
    dadx *= 1. / L2Norm(dadx);
    return dadx;
}

Eigen::VectorXd ddcDSI::zdiff(const Eigen::VectorXd& a) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    vector2field(a, u, temp, salt);
    Eigen::VectorXd dadz(a.size());
    dadz.setZero();
    u = chflow::zdiff(u);
    temp = chflow::zdiff(temp);
    salt = chflow::zdiff(salt);
    field2vector(u, temp, salt, dadz);
    dadz *= 1. / L2Norm(dadz);
    return dadz;
}

Eigen::VectorXd ddcDSI::tdiff(const Eigen::VectorXd& a, Real epsDt) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField temp(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField salt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FieldSymmetry sigma;
    Real T;
    // quick hack to avoid new interface or creating simple f() for DDC
    TimeStep dt = TimeStep(epsDt, 0, 1, 1, 0, 0, false);
    int fcount = 0;
    PoincareCondition* h = 0;
    Real CFL = 0.0;
    FlowField edudtf(u);
    FlowField edtempdtf(temp);
    FlowField edsaltdtf(salt);
    std::ostream muted_os(0);
    //   vector2field (a, u, temp);
    extractVectorDDC(a, u, temp, salt, sigma, T);
    // use existing f() instead of simple
    f(u, temp, salt, epsDt, h, edudtf, edtempdtf, edsaltdtf, ddcflags_, dt, fcount, CFL, muted_os);
    //   f (temp, 1,epsDt, edtempdtf, dnsflags_, *os_);
    edudtf -= u;
    edtempdtf -= temp;
    edsaltdtf -= salt;
    Eigen::VectorXd dadt(a.size());
    field2vector(edudtf, edtempdtf, edsaltdtf, dadt);
    dadt *= 1. / L2Norm(dadt);
    return dadt;
}

void ddcDSI::updateMu(Real mu) {
    DSI::updateMu(mu);
    if (ddc_cPar_ == ddc_continuationParameter::none) {
        cfDSI::updateMu(mu);
    } else if (ddc_cPar_ == ddc_continuationParameter::kxkz) {
        ddcflags_.kxkz = mu;
        Lx_ = 2*M_PI/mu;
        Lz_ = 2*M_PI/mu;
    } else if (ddc_cPar_ == ddc_continuationParameter::Lx) {
        Lx_ = mu;
    } else if (ddc_cPar_ == ddc_continuationParameter::Lz) {
        Lz_ = mu;
    }else if (ddc_cPar_ == ddc_continuationParameter::Ri) {
        ddcflags_.Ri = mu;
    }else if (ddc_cPar_ == ddc_continuationParameter::Rrho) {
        ddcflags_.Rrho = mu;
    }else if (ddc_cPar_ == ddc_continuationParameter::Pr) {
        ddcflags_.Pr = mu;
    }else if (ddc_cPar_ == ddc_continuationParameter::Ra) {
        ddcflags_.Ra = mu;
    }else if (ddc_cPar_ == ddc_continuationParameter::Le) {
        ddcflags_.Le = mu;
    }else if (ddc_cPar_ == ddc_continuationParameter::Rsep) {
        ddcflags_.Rsep = mu;
    }else if (ddc_cPar_ == ddc_continuationParameter::Rey) {
        ddcflags_.Rey = mu;
    }else if (ddc_cPar_ == ddc_continuationParameter::gammax) {
        ddcflags_.gammax = mu;
    }else if (ddc_cPar_ == ddc_continuationParameter::gammaz) {
        ddcflags_.gammaz = mu;
    }else if (ddc_cPar_ == ddc_continuationParameter::Uw) {
        ddcflags_.Uwall = mu;
        ddcflags_.ulowerwall = -0.5*mu;
        ddcflags_.uupperwall = 0.5*mu;
    }else if (ddc_cPar_ == ddc_continuationParameter::Ek) {
        ddcflags_.Ek = mu;
    }else {
        throw invalid_argument("ddcDSI::updateMu(): continuation parameter is unknown");
    }
}

void ddcDSI::chooseMuDDC(string muName) {
    ddc_continuationParameter ddc_cPar = s2ddc_cPar(muName);

    if (ddc_cPar == ddc_continuationParameter::none)
        ddcDSI::chooseMu(muName);
    else
        chooseMuDDC(ddc_cPar);
}

void ddcDSI::chooseMuDDC(ddc_continuationParameter mu) {
    ddc_cPar_ = mu;
    switch (mu) {
        case ddc_continuationParameter::kxkz:
            updateMu(ddcflags_.kxkz);
            break;
        case ddc_continuationParameter::Lx:
            updateMu(Lx_);
            break;
        case ddc_continuationParameter::Lz:
            updateMu(Lz_);
            break;
        case ddc_continuationParameter::Ri:
            updateMu(ddcflags_.Ri);
            break;
        case ddc_continuationParameter::Rrho:
            updateMu(ddcflags_.Rrho);
            break;
        case ddc_continuationParameter::Pr:
            updateMu(ddcflags_.Pr);
            break;
        case ddc_continuationParameter::Ra:
            updateMu(ddcflags_.Ra);
            break;
        case ddc_continuationParameter::Le:
            updateMu(ddcflags_.Le);
            break;
        case ddc_continuationParameter::Rsep:
            updateMu(ddcflags_.Rsep);
            break;
        case ddc_continuationParameter::Rey:
            updateMu(ddcflags_.Rey);
            break;
        case ddc_continuationParameter::gammax:
            updateMu(ddcflags_.gammax);
            break;
        case ddc_continuationParameter::gammaz:
            updateMu(ddcflags_.gammaz);
            break;
        case ddc_continuationParameter::Uw:
            updateMu(ddcflags_.Uwall);
            break;
        case ddc_continuationParameter::Ek:
            updateMu(ddcflags_.Ek);
            break;
        case ddc_continuationParameter::none:
            throw invalid_argument(
                "ddcDSI::chooseMu(): continuation parameter is none, we should not reach this point");
        default:
            throw invalid_argument("ddcDSI::chooseMu(): continuation parameter is unknown");
    }
}

ddc_continuationParameter ddcDSI::s2ddc_cPar(string muname) {
    std::transform(muname.begin(), muname.end(), muname.begin(), ::tolower);  // why is the string made lower case?
    if (muname == "kxkz")
        return ddc_continuationParameter::kxkz;
    else if (muname == "lx")
        return ddc_continuationParameter::Lx;
    else if (muname == "lz")
        return ddc_continuationParameter::Lz;
    else if (muname == "ri")
        return ddc_continuationParameter::Ri;
    else if (muname == "rrho")
        return ddc_continuationParameter::Rrho;
    else if (muname == "pr")
        return ddc_continuationParameter::Pr;
    else if (muname == "ra")
        return ddc_continuationParameter::Ra;
    else if (muname == "le")
        return ddc_continuationParameter::Le;
    else if (muname == "rsep")
        return ddc_continuationParameter::Rsep;
    else if (muname == "rey")
        return ddc_continuationParameter::Rey;
    else if (muname == "gammax")
        return ddc_continuationParameter::gammax;
    else if (muname == "gammaz")
        return ddc_continuationParameter::gammaz;
    else if (muname == "uw")
        return ddc_continuationParameter::Uw;
    else if (muname == "ek")
        return ddc_continuationParameter::Ek;
    else {
        cout << "ddcDSI::s2ddc_cPar(): ddc_continuation parameter '"+muname+"' is unknown, defaults to 'none'" << endl;
        return ddc_continuationParameter::none;
    }
}

string ddcDSI::printMu() { return ddc_cPar2s(ddc_cPar_); }

string ddcDSI::ddc_cPar2s(ddc_continuationParameter ddc_cPar) {
    if (ddc_cPar == ddc_continuationParameter::none)
        return cfDSI::cPar2s(cPar_);
    else if (ddc_cPar == ddc_continuationParameter::kxkz)
        return "kxkz";
    else if (ddc_cPar == ddc_continuationParameter::Lx)
        return "Lx";
    else if (ddc_cPar == ddc_continuationParameter::Lz)
        return "Lz";
    else if (ddc_cPar == ddc_continuationParameter::Ri)
        return "Ri";
    else if (ddc_cPar == ddc_continuationParameter::Rrho)
        return "Rrho";
    else if (ddc_cPar == ddc_continuationParameter::Pr)
        return "Pr";
    else if (ddc_cPar == ddc_continuationParameter::Ra)
        return "Ra";
    else if (ddc_cPar == ddc_continuationParameter::Le)
        return "Le";
    else if (ddc_cPar == ddc_continuationParameter::Rsep)
        return "Rsep";
    else if (ddc_cPar == ddc_continuationParameter::Rey)
        return "Rey";
    else if (ddc_cPar == ddc_continuationParameter::gammax)
        return "gammax";
    else if (ddc_cPar == ddc_continuationParameter::gammaz)
        return "gammaz";
    else if (ddc_cPar == ddc_continuationParameter::Uw)
        return "Uw";
    else if (ddc_cPar == ddc_continuationParameter::Ek)
        return "Ek";
    else
        throw invalid_argument("ddcDSI::ddc_cPar2s(): continuation parameter is not convertible to string");
}

void ddcDSI::saveParameters(string searchdir) {
    // cfDSI::saveParameters (searchdir);
    ddcflags_.save(searchdir);
}

void ddcDSI::saveEigenvec(const Eigen::VectorXd& ev, const string label, const string outdir) {
    FlowField efu(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField eft(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField efs(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    vector2field(ev, efu, eft, efs);
    efu *= 1.0 / L2Norm(efu);
    eft *= 1.0 / L2Norm(eft);
    efs *= 1.0 / L2Norm(efs);
    efu.save(outdir + "efu" + label);
    eft.save(outdir + "eft" + label);
    efs.save(outdir + "efs" + label);
}

void ddcDSI::saveEigenvec(const Eigen::VectorXd& evA, const Eigen::VectorXd& evB, const string label1,
                          const string label2, const string outdir) {
    FlowField efAu(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField efBu(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField efAt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField efBt(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField efAs(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField efBs(Nx_, Ny_, Nz_, 1, Lx_, Lz_, ya_, yb_, cfmpi_);
    vector2field(evA, efAu, efAt, efAs);
    vector2field(evB, efBu, efBt, efBs);
    Real cu = 1.0 / sqrt(L2Norm2(efAu) + L2Norm2(efBu));
    Real ct = 1.0 / sqrt(L2Norm2(efAt) + L2Norm2(efBt));
    Real cs = 1.0 / sqrt(L2Norm2(efAs) + L2Norm2(efBs));
    efAu *= cu;
    efBu *= cu;
    efAt *= ct;
    efBt *= ct;
    efAs *= cs;
    efBs *= cs;
    efAu.save(outdir + "efu" + label1);
    efBu.save(outdir + "efu" + label2);
    efAt.save(outdir + "eft" + label1);
    efBt.save(outdir + "eft" + label2);
    efAs.save(outdir + "efs" + label1);
    efBs.save(outdir + "efs" + label2);
}

/* OUTSIDE CLASS */

// G(x) = G(u,sigma) = (sigma f^T(u) - u) for orbits
void G(const FlowField& u, const FlowField& temp, const FlowField& salt, Real& T, PoincareCondition* h, const FieldSymmetry& sigma,
       FlowField& Gu, FlowField& Gtemp, FlowField& Gsalt, const DDCFlags& ddcflags, const TimeStep& dt, bool Tnormalize, Real Unormalize,
       int& fcount, Real& CFL, ostream& os) {
    f(u, temp, salt, T, h, Gu, Gtemp, Gsalt, ddcflags, dt, fcount, CFL, os);
    Real funorm = L2Norm3d(Gu);
    Gu *= sigma;
    Gu -= u;
    if (sigma.sy() == -1) {
        // wall-normal mirroring in velocity requires sign change in temperature and salinity
        FieldSymmetry tsigma(sigma);
        FieldSymmetry ssigma(sigma);
        FieldSymmetry inv(-1);
        tsigma *= inv;
        ssigma *= inv;
        Gtemp *= tsigma;
        Gsalt *= ssigma;
    } else {
        Gtemp *= sigma;
        Gsalt *= sigma;
    }
    Gtemp -= temp;
    Gsalt -= salt;

    if (Tnormalize) {
        Gu *= 1.0 / T;
        Gtemp *= 1.0 / T;
        Gsalt *= 1.0 / T;
    }
    if (Unormalize != 0.0) {
        Gu *= 1. / sqrt(abs(funorm * (Unormalize - funorm)));
        // u should stay off zero, so normalize with u for now - temp should also stay away from zero
        Gtemp *= 1. / sqrt(abs(funorm * (Unormalize - funorm)));
        Gsalt *= 1. / sqrt(abs(funorm * (Unormalize - funorm)));
    }
}

void f(const FlowField& u, const FlowField& temp, const FlowField& salt, Real& T, PoincareCondition* h, FlowField& f_u, FlowField& f_temp, FlowField& f_salt,
       const DDCFlags& ddcflags_, const TimeStep& dt_, int& fcount, Real& CFL, ostream& os) {
    if (!isfinite(L2Norm(u))) {
        os << "error in f: u is not finite. exiting." << endl;
        exit(1);
    }
    DDCFlags flags(ddcflags_);
    flags.logstream = &os;
    TimeStep dt(dt_);
    vector<FlowField> fields = {u, temp, salt, FlowField(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi())};

    //   f_u = u;
    //   f_temp = temp;
    //   f_salt = salt;
    // No Poincare section, just integration to time T
    if (h == 0) {
        if (T < 0) {
            os << "f: negative integration time T == " << T << endl
               << "returning f(u,T) == (1+abs(T))*u" << endl
               << "returning f(temp,T) == (1+abs(T))*temp" << endl
               << "returning f(salt,T) == (1+abs(T))*salt" << endl;
            fields[0] *= 1 + abs(T);
            fields[1] *= 1 + abs(T);
            fields[2] *= 1 + abs(T);
            return;
        }
        // Special case #1: no time steps
        if (T == 0) {
            os << "f: T==0, no integration, returning u, temp, and salt" << endl;
            return;
        }
        dt.adjust_for_T(T, false);
        flags.dt = dt;
        // Adjust dt for CFL if necessary
        DDC ddc(fields, flags);
        ddc.advance(fields, 1);
        if (dt.variable()) {
            dt.adjust(ddc.CFL(fields[0]), false);
            ddc.reset_dt(dt);
        }
        //  t == current time in integration
        //  T == total integration time
        // dT == CFL-check interval
        // dt == DNS time-step
        //  N == T/dT,  n == dT/dt;
        //  T == N dT, dT == n dt
        //  t == s dT (s is loop index)

        os << "f^T: " << flush;
        for (int s = 1; s <= dt.N(); ++s) {
            Real t = s * dt.dT();
            CFL = ddc.CFL(fields[0]);
            if (s % 10 == 0)
                os << iround(t) << flush;
            else if (s % 2 == 0) {
                if (CFL < dt.CFLmin())
                    os << '<' << flush;
                else if (CFL > dt.CFLmax())
                    os << '>' << flush;
                else
                    os << '.' << flush;
            }
            ddc.advance(fields, dt.n());
            if (dt.variable() && dt.adjust(CFL, false))
                ddc.reset_dt(dt);
        }

    }
    // Poincare section computation: return Poincare crossing nearest to t=T, with Tmin < t < Tmax.
    else {
        cout << "Poincare sectioning not yet implemented (markd as experimental)." << endl;
        exit(1);
        /*    // Adjust dt for CFL if necessary
            DNSPoincare dns (f_u, h, flags);
            if (dt.variable()) {
              dns.advance (f_u, p, 1);
              dt.adjust (dns.CFL());
              dns.reset_dt (dt,u);
              f_u = u;
            }
            // Collect all Poincare crossings between Tfudgemin and Tfudgemax
            // If we don't find one in that range, go as far as Tlastchance
            Real dTfudge = 1.05;
            Real Tfudgemin = lesser (0.90*T, T - dTfudge*dt.dT());
            Real Tfudgemax = Greater (1.02*T, T + dTfudge*dt.dT());
            Real Tlastchance = 10*T;

            vector<FlowField> ucross;
            vector<Real> tcross;
            vector<int>  scross;
            int s=0;
            int crosssign = 0; // look for crossings in either direction

            os << "f^t: " << flush;

            for (Real t=0; t<=Tlastchance; t += dt.dT(), ++s) {

              CFL = dns.CFL();

              if (s % 10 == 0)   os << iround (t);
              else if (s % 2 == 0) {
                if (CFL > dt.CFLmax())  os << '>';
                else if (CFL < dt.CFLmin())  os << '<';
                else  os << '.';
                os << flush;
              }

              // Collect any Poincare crossings
              bool crossed = dns.advanceToSection (f_u, p, dt.n(), crosssign, Tfudgemin);
              if (crossed && t >= Tfudgemin) {
                ucross.push_back (dns.ucrossing());
                tcross.push_back (dns.tcrossing());
                scross.push_back (dns.scrossing());
              }

              // If we've found at least one crossing within the fudge range, stop.
              // Otherwise continue trying until Tlastchance
              if (ucross.size() > 0 && t >= Tfudgemax)
                break;
            }

            if (ucross.size() <1) {

              os << "\nError in f(u, T, f_u, flags, dt, fcount, CFL, os) :\n";
              os << "the integration did not reach the Poincare section.\n";
              os << "Returning laminar solution and a b.s. value for the crossing time.\n";
              os << "I hope you can do something useful with them." << endl;
              f_u.setToZero();
              T = dns.time();
              ++fcount;
              return;
            }
            os << "  " << flush;

            // Now select the crossing that is closest to the estimated crossing time
            FlowField ubest = ucross[0];
            Real  Tbest = tcross[0];
            int   sbest = scross[0];
            int   nbest = 0;

            for (uint n=1; n<ucross.size(); ++n) {
              if (abs (tcross[n]-T) < abs (Tbest-T)) {
                ubest = ucross[n];
                Tbest = tcross[n];
                sbest = scross[n];
                nbest = n;
              }
            }
            os << nbest << (sbest==1 ? '+' : '-') << " at t== " << Tbest << flush;

            T = Tbest;
            f_u = ubest;

            // Now check if there are any crossings of opposite sign close by.
            // This signals near-tangency to Poincare section, which'll mess up
            // the search. Just print warning and let user intervene manually.
            Real l2distubestucross = 0;
            for (uint n=0; n<ucross.size(); ++n) {
              l2distubestucross = L2Dist (ubest, ucross[n]);
              if ( (u.taskid() == 0) && (scross[n] != sbest)) {
                os << "\nWARNING : There is a nearby Poincare crossing of opposite sign," << endl;
                os << "signalling near-tangency to section. You should probably switch " << endl;
                os << "to another Poincare crossing." << endl;
                os << "(ubest, unear) signs == " << sbest << ", " << scross[n] << endl;
                os << "(ubest, unear) times == " << Tbest << ", " << tcross[n] << endl;
                os << "(ubest, unear) dist  == " << l2distubestucross << endl;
              }
            }*/
    }

    if (!isfinite(L2Norm(f_u))) {
        os << "error in f: f(u,t) is not finite. exiting." << endl;
        exit(1);
    }
    if (!isfinite(L2Norm(f_temp))) {
        os << "error in f: f(temp,t) is not finite. exiting." << endl;
        exit(1);
    }
    if (!isfinite(L2Norm(f_salt))) {
        os << "error in f: f(salt,t) is not finite. exiting." << endl;
        exit(1);
    }

    ++fcount;
    f_u = fields[0];
    f_temp = fields[1];
    f_salt = fields[2];
    return;
}

}  // namespace chflow