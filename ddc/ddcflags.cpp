/**
 *
 * Original author: Duc Nguyen
 */

#include "modules/ddc/macros.h"
#include "modules/ddc/ddcflags.h"
using namespace std;
namespace chflow {

DDCFlags::DDCFlags(Real kxkz_,
    Real Rey_, Real Pr_, Real Ra_, Real Le_, Real Rrho_, Real Rsep_, Real Ri_, Real Ek_, Real gammax_, Real gammaz_,
                   Real ulowerwall_, Real uupperwall_, 
                   Real wlowerwall_, Real wupperwall_,
                   Real tlowerwall_, Real tupperwall_, 
                   Real slowerwall_, Real supperwall_, 
                   Real ystats_)
    : kxkz(kxkz_),
      Rey(Rey_),
      Pr(Pr_),
      Ra(Ra_),
      Le(Le_),
      Rrho(Rrho_),
      Rsep(Rsep_),
      Ri(Ri_),  
      Ek(Ek_),  
      gammax(gammax_),
      gammaz(gammaz_),
      
      
      tlowerwall(tlowerwall_),
      tupperwall(tupperwall_),
      slowerwall(slowerwall_),
      supperwall(supperwall_),

      
      ystats(ystats_) {
    
    ulowerwall = ulowerwall_;
    uupperwall = uupperwall_;
    wlowerwall = wlowerwall_;
    wupperwall = wupperwall_;

    // timestepping = CNRK2;
    // constraint = PressureGradient;
    // baseflow = LaminarBase;
    // nonlinearity = Convection;
    // dealiasing = DealiasXZ;
}

DDCFlags::DDCFlags(ArgList& args, const bool laurette) {
    // DDC system parameters
    args.section("System parameters");
    const Real kxkz_ = args.getreal("-kxkz", "--kxkz", 1, "wavenumber in x and z directions");
    const Real Rey_ = args.getreal("-Rey", "--Reynolds", 400, "pseudo-Reynolds number == 1/nu");
    const Real Pr_ = args.getreal("-Pr", "--Prandtl", 10, "Prandtl number == nu/kpT");
    const Real Ra_ = args.getreal("-Ra", "--Rayleigh", 1000, "Thermal Rayleigh number");
    const Real Le_ = args.getreal("-Le", "--Lewis", 100, "Lewis number");
    const Real Rrho_ = args.getreal("-Rr", "--Rrho", 2, "Density stability ratio");
    const Real Rsep_ = args.getreal("-Rs", "--Rsep", 1, "Separation ratio for binary fluid convection");
    const Real Ri_ = args.getreal("-Ri", "--Richardson", 10, "Richardson number");
    const Real Ek_ = args.getreal("-Ek", "--Ekman", 1e-1, "Ekman number");
    const Real gammax_ = args.getreal("-GammaX", "--GammaX", 0, "Gamma angle in X");
    const Real gammaz_ = args.getreal("-GammaZ", "--GammaZ", 0, "Gamma angle in Z");
    
    // define Channelflow boundary conditions from arglist
    args2BC(args);

    // define Channelflow numerics from arglist
    args2numerics(args, laurette);

    const std::string tsymmstr = args.getstr("-tsymms", "--tempsymmetries", "",
                                             "constrain temp(t) to invariant "
                                             "symmetric subspace, argument is the filename for a file "
                                             "listing the generators of the isotropy group");
    const std::string ssymmstr = args.getstr("-ssymms", "--saltsymmetries", "",
                                             "constrain salt(t) to invariant "
                                             "symmetric subspace, argument is the filename for a file "
                                             "listing the generators of the isotropy group");
    const Real ystats_ = args.getreal("-ys", "--ystats", 0, "y-coordinate of height dependent statistics, e.g. Nu(y)");
    
    // set flags
    ystats = ystats_;
    kxkz = kxkz_;
    Rey = Rey_;
    Pr = Pr_;
    Ra = Ra_;
    Le = Le_;
    Rrho = Rrho_;
    Rsep = Rsep_;
    Ri = Ri_;
    Ek = Ek_; // Ekman number, used for rotating RBC
    gammax = gammax_ / 180.0 * pi;
    gammaz = gammaz_ / 180.0 * pi;
    
    nu = P1; // get kinematic viscosity of system


    const Real ulowerwall_ = args.getreal("-Ua", "--ulowerwall", 0, "X-Velocity at lower wall, U(y=a)");
    const Real uupperwall_ = args.getreal("-Ub", "--uupperwall", 0, "X-Velocity at upper wall, U(y=b)");
    const Real wlowerwall_ = args.getreal("-Wa", "--wlowerwall", 0, "Z-Velocity at lower wall, W(y=a)");
    const Real wupperwall_ = args.getreal("-Wb", "--wupperwall", 0, "Z-Velocity at upper wall, W(y=b)");

    ulowerwall = ulowerwall_;
    uupperwall = uupperwall_;
    wlowerwall = wlowerwall_;
    wupperwall = wupperwall_;

    const Real tlowerwall_ = args.getreal("-Ta", "--tlowerwall", 0, "Temperature at lower wall, T(y=a)");
    const Real tupperwall_ = args.getreal("-Tb", "--tupperwall", 1, "Temperature at upper wall, T(y=b)");
    const Real slowerwall_ = args.getreal("-Sa", "--slowerwall", 0, "Salinity at lower wall, S(y=a)");
    const Real supperwall_ = args.getreal("-Sb", "--supperwall", 1, "Salinity at upper wall, S(y=b)");

    tlowerwall = tlowerwall_;
    tupperwall = tupperwall_;
    slowerwall = slowerwall_;
    supperwall = supperwall_;

    // timestepping = CNRK2;
    // constraint = PressureGradient;
    // baseflow = LaminarBase;
    // nonlinearity = Convection;
    // dealiasing = DealiasXZ;
    
    if (tsymmstr.length() > 0) {
        SymmetryList tsymms(tsymmstr);
        tempsymmetries = tsymms;
    }
    if (ssymmstr.length() > 0) {
        SymmetryList ssymms(ssymmstr);
        saltsymmetries = ssymms;
    }
    save();
}

void DDCFlags::save(const std::string& savedir) const {
    DNSFlags::save(savedir);
    if (mpirank() == 0) {
        std::string filename = appendSuffix(savedir, "ddcflags.txt");
        std::ofstream os(filename.c_str());
        if (!os.good())
            cferror("DDCFlags::save(savedir) :  can't open file " + filename);
        os.precision(16);
        os.setf(std::ios::left);
        os << std::setw(REAL_IOWIDTH) << Rey << "  %Rey\n"
           << std::setw(REAL_IOWIDTH) << Pr << "  %Pr\n"
           << std::setw(REAL_IOWIDTH) << Ra << "  %Ra\n"
           << std::setw(REAL_IOWIDTH) << Le << "  %Le\n"
           << std::setw(REAL_IOWIDTH) << Rrho << "  %Rrho\n"
           << std::setw(REAL_IOWIDTH) << Rsep << "  %Rsep\n"
           << std::setw(REAL_IOWIDTH) << Ri << "  %Ri\n"
           << std::setw(REAL_IOWIDTH) << Ek << "  %Ek\n"
           << std::setw(REAL_IOWIDTH) << gammax << "  %GammaX\n"
           << std::setw(REAL_IOWIDTH) << gammaz << "  %GammaZ\n"
           << std::setw(REAL_IOWIDTH) << uupperwall << "  %uupperwall\n"
           << std::setw(REAL_IOWIDTH) << ulowerwall << "  %ulowerwall\n"
           << std::setw(REAL_IOWIDTH) << wupperwall << "  %wupperwall\n"
           << std::setw(REAL_IOWIDTH) << wlowerwall << "  %wlowerwall\n"
           << std::setw(REAL_IOWIDTH) << tupperwall << "  %tupperwall\n"
           << std::setw(REAL_IOWIDTH) << tlowerwall << "  %tlowerwall\n"
           << std::setw(REAL_IOWIDTH) << supperwall << "  %supperwall\n"
           << std::setw(REAL_IOWIDTH) << slowerwall << "  %slowerwall\n"
           << std::setw(REAL_IOWIDTH) << ystats << "  %ystats\n";
        os.unsetf(std::ios::left);
    }
}

void DDCFlags::load(int taskid, const std::string indir) {
    DNSFlags::load(taskid, indir);
    std::ifstream is;
    if (taskid == 0) {
        is.open(indir + "ddcflags.txt");
        if (!is.good())
            cferror(" DDCFlags::load(taskid, flags, dt, indir):  can't open file " + indir + "ddcflags.txt");
    }
    Rey = getRealfromLine(taskid, is);
    Pr = getRealfromLine(taskid, is);
    Ra = getRealfromLine(taskid, is);
    Le = getRealfromLine(taskid, is);
    Rrho = getRealfromLine(taskid, is);
    Rsep = getRealfromLine(taskid, is);
    Ri = getRealfromLine(taskid, is);
    Ek = getRealfromLine(taskid, is);
    gammax = getRealfromLine(taskid, is);
    gammaz = getRealfromLine(taskid, is);
    uupperwall = getRealfromLine(taskid, is);
    ulowerwall = getRealfromLine(taskid, is);
    wupperwall = getRealfromLine(taskid, is);
    wlowerwall = getRealfromLine(taskid, is);
    tupperwall = getRealfromLine(taskid, is);
    tlowerwall = getRealfromLine(taskid, is);
    supperwall = getRealfromLine(taskid, is);
    slowerwall = getRealfromLine(taskid, is);
    ystats = getRealfromLine(taskid, is);
}

}  // namespace chflow
