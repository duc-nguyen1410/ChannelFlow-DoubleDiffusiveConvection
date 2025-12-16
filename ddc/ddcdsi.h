/**
 * Dynamical System Interface for the DDC module
 *
 * Original author: Duc Nguyen
 */

#ifndef DDCDSI_H
#define DDCDSI_H

#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "cfbasics/cfvector.h"
#include "channelflow/cfdsi.h"
#include "channelflow/cfmpi.h"
#include "channelflow/chebyshev.h"
#include "channelflow/flowfield.h"
#include "channelflow/periodicfunc.h"
#include "channelflow/symmetry.h"
#include "channelflow/tausolver.h"
#include "channelflow/utilfuncs.h"
#include "modules/ddc/ddc.h"
#include "nsolver/nsolver.h"

using namespace std;

namespace chflow {

enum class ddc_continuationParameter {
    none,
    kxkz,
    Lx,
    Lz,
    Ri,
    Rrho,
    Pr,
    Ra,
    Le,
    Rsep,
    Rey,
    gammax,
    gammaz,
    Uw,
    Ek
};

// Real GMRESHookstep_vector (FlowField& u, FlowField& alpha, Real& T, FieldSymmetry& sigma,
//                            PoincareCondition* h,
//                            const nsolver::hookstepSearchFlags& searchflags,
//                            DNSFlags& dnsflags, VEDNSFlags& vednsflags, TimeStep& dt, Real& CFL, Real Unormalize);

std::vector<Real> ddcstats(const FlowField& u, const FlowField& temp, const FlowField& salt, const DDCFlags flags = DDCFlags());

// header for fieldstats
string ddcfieldstatsheader(const DDCFlags flags = DDCFlags());

// header for fieldstats with parameter t
string ddcfieldstatsheader_t(std::string muname, const DDCFlags flags = DDCFlags());
string ddcfieldstats(const FlowField& u, const FlowField& temp, const FlowField& salt, const DDCFlags flags = DDCFlags());
string ddcfieldstats_t(const FlowField& u, const FlowField& temp, const FlowField& salt, Real t, const DDCFlags flags = DDCFlags());
FlowField totalVelocity(const FlowField& velo, const DDCFlags flags);
FlowField totalTemperature(const FlowField& temp, const DDCFlags flags);
FlowField totalSalinity(const FlowField& salt, const DDCFlags flags);
Real heatcontent(const FlowField& ttot, const DDCFlags flags);
Real saltcontent(const FlowField& stot, const DDCFlags flags);
Real dissipation(const FlowField& u, const DDCFlags flags, bool normalize = true, bool relative = true);
Real heatflux(const FlowField& ttot, const DDCFlags flags, bool normalize = true, bool relative = true);
Real massflux(const FlowField& stot, const DDCFlags flags, bool normalize = true, bool relative = true);
Real Nusselt_t_plane(const FlowField& utot, const FlowField& ttot, const DDCFlags flags, bool relative = true);
Real Nusselt_c_plane(const FlowField& utot, const FlowField& stot, const DDCFlags flags, bool relative = true);
Real buoyPowerInput(const FlowField& utot, const FlowField& ttot, const FlowField& stot, const DDCFlags flags, bool relative=true);
// Real UdPcontent(const FlowField& utot, const DDCFlags flags);
class ddcDSI : public cfDSI {
   public:
    /** \brief default constructor */
    ddcDSI();
    virtual ~ddcDSI() {}

    /** \brief Initialize ddcDSI */
    ddcDSI(DDCFlags& ddcflags, FieldSymmetry sigma, PoincareCondition* h, TimeStep dt, bool Tsearch, bool xrelative,
           bool zrelative, bool Tnormalize, Real Unormalize, const FlowField& u, const FlowField& temp, const FlowField& salt,
           std::ostream* os = &std::cout);

    Eigen::VectorXd eval(const Eigen::VectorXd& x) override;
    Eigen::VectorXd eval(const Eigen::VectorXd& x0, const Eigen::VectorXd& x1,
                         bool symopt) override;  // needed for multishooting
    void save(const Eigen::VectorXd& x, const string filebase, const string outdir = "./",
              const bool fieldsonly = false) override;

    string stats(const Eigen::VectorXd& x) override;
    pair<string, string> stats_minmax(const Eigen::VectorXd& x) override;
    string statsHeader() override;
    void phaseShift(Eigen::VectorXd& x) override;
    void phaseShift(Eigen::MatrixXd& y) override;
    Real extractT(const Eigen::VectorXd& x) override;
    Real extractXshift(const Eigen::VectorXd& x) override;
    Real extractZshift(const Eigen::VectorXd& x) override;

    void makeVectorDDC(const FlowField& u, const FlowField& temp, const FlowField& salt, const FieldSymmetry& sigma, const Real T,
                       Eigen::VectorXd& x);
    void extractVectorDDC(const Eigen::VectorXd& x, FlowField& u, FlowField& temp, FlowField& salt, FieldSymmetry& sigma, Real& T);

    /// \name Compute derivatives of the two FlowFields contained in this vector
    Eigen::VectorXd xdiff(const Eigen::VectorXd& a) override;
    Eigen::VectorXd zdiff(const Eigen::VectorXd& a) override;
    Eigen::VectorXd tdiff(const Eigen::VectorXd& a, Real epsDt) override;
    /// \name Handle continuation parameter
    void updateMu(Real mu) override;
    void chooseMuDDC(std::string muName);
    void chooseMuDDC(ddc_continuationParameter mu);
    string printMu() override;  // document
    void saveParameters(string searchdir) override;
    ddc_continuationParameter s2ddc_cPar(std::string muName);
    string ddc_cPar2s(ddc_continuationParameter cPar);

    // Save real eigenvectors
    void saveEigenvec(const Eigen::VectorXd& x, const string label, const string outdir) override;
    // Save complex conjugate eigenvectors pair
    void saveEigenvec(const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const string label1, const string label2,
                      const string outdir) override;

   protected:
    DDCFlags ddcflags_;
    ddc_continuationParameter ddc_cPar_ = ddc_continuationParameter::none;
};

// G(x) = G(u,sigma) = (sigma f^T(u) - u) for orbits
void G(const FlowField& u, const FlowField& temp, const FlowField& salt, Real& T, PoincareCondition* h, const FieldSymmetry& sigma,
       FlowField& Gu, FlowField& Gtemp, FlowField& Gsalt, const DDCFlags& ddcflags, const TimeStep& dt, bool Tnormalize, Real Unormalize,
       int& fcount, Real& CFL, ostream& os);
void f(const FlowField& u, const FlowField& temp, const FlowField& salt, Real& T, PoincareCondition* h, FlowField& f_u, FlowField& f_temp, FlowField& f_salt,
       const DDCFlags& ddcflags_, const TimeStep& dt_, int& fcount, Real& CFL, ostream& os);

}  // namespace chflow

#endif
