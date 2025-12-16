/**
 * Control parameters for time-integration within the DDC module
 * DDCFlags specifies all relevant parameters for integrating the double-diffusive equations in doubly periodic
 * channel domains.
 *
 * Original author: Duc Nguyen
 */

#ifndef DDCFLAGS_H
#define DDCFLAGS_H

#include "channelflow/dnsflags.h"
#include "channelflow/utilfuncs.h"
#include "modules/ddc/ddcflags.h"
namespace chflow {

/** \brief extension of the DNSFlags class for DDC
 *
 * DDCFlags class, holds all additional parameters for convective shear flows
 */
class DDCFlags : public DNSFlags {
    // Is derived from DNSFlags to keep saving and loading consistent

   public:
    DDCFlags(Real kxkz = 1.0,
        Real Rey = 400.0, Real Pr = 10.0, Real Ra = 1000.0, Real Le = 100.0, Real Rrho = 2.0, Real Rsep = 1.0, Real Ri = 10.0, Real Ek=1e-5, Real gammax = 0.0, Real gammaz = 0.0, 
             Real ulowerwall = 0.0, Real uupperwall = 0.0, 
             Real wlowerwall = 0.0, Real wupperwall = 0.0, 
             Real tlowerwall = 0.0, Real tupperwall = 1.0, 
             Real slowerwall = 0.0, Real supperwall = 1.0, 
             Real ystats = 0);
    //     DDCFlags (const DDCFlags& ddcflags);
    //     DDCFlags (const DNSFlags& flags);
    //     DDCFlags (const std::string& filebase);
    DDCFlags(ArgList& args, const bool laurette = false);

    /** \brief The infamous virtual destructor */
    virtual ~DDCFlags() = default;

    Real kxkz;

    Real Rey;
    Real Pr;
    Real Ra;
    Real Le;
    Real Rrho;
    Real Rsep; // use only for binary fluid convection
    Real Ri;
    Real Ek; // Ekman number, used for rotating flows
    Real gammax;
    Real gammaz;

    Real tlowerwall;
    Real tupperwall;
    Real slowerwall;
    Real supperwall;

    Real ystats;

    cfarray<FieldSymmetry> tempsymmetries;  // restrict temp(t) to these symmetries
    cfarray<FieldSymmetry> saltsymmetries;
    
    // override DNSFlags::save/load. Both methods include a call to the parent method
    virtual void save(const std::string& savedir = "") const override;
    virtual void load(int taskid, const std::string indir) override;
};

}  // namespace chflow
#endif  // DDCFLAGS_H