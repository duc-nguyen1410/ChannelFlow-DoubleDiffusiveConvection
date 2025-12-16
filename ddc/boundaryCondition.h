/**
 * free-slip boundary condition on walls
 *
 * Original author: Duc Nguyen
 */

#ifndef FREESLIPBC_H
#define FREESLIPBC_H

#include "channelflow/diffops.h"
#include "channelflow/flowfield.h"
#include "modules/ddc/macros.h"
#ifdef FREESLIP
namespace chflow {
void freeslipBC(std::vector<FlowField>& fields) {
    
    if (fields[0].Nd() != 3) {
        cferror("freeslipBC(u): u.Nd() = " + i2s(fields[0].Nd()) + " != 3");
        exit(1);
    }
    
    fields[0].makePhysical();

    lint Nz = fields[0].Nz();
    lint nxlocmin = fields[0].nxlocmin();
    lint nxlocmax = fields[0].nxlocmin() + fields[0].Nxloc();
    lint nylocmin = fields[0].nylocmin();
    lint nylocmax = fields[0].nylocmax();
    
    // top
    if(nylocmin == 0){
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx){
            for (lint nz = 0; nz < Nz; ++nz) {
                fields[0](nx, nylocmin, nz, 0) = fields[0](nx, nylocmin+1, nz, 0); // du/dy=0
                fields[0](nx, nylocmin, nz, 1) = 0; // No-penetration
                fields[0](nx, nylocmin, nz, 2) = fields[0](nx, nylocmin+1, nz, 2); // dw/dy=0
            }
        }
    }

    // bottom
    if(nylocmax == fields[0].Ny()-1){
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx){
            for (lint nz = 0; nz < Nz; ++nz) {
                fields[0](nx, nylocmax, nz, 0) = fields[0](nx, nylocmax-1, nz, 0); // du/dy=0
                fields[0](nx, nylocmax, nz, 1) = 0; // No-penetration
                fields[0](nx, nylocmax, nz, 2) = fields[0](nx, nylocmax-1, nz, 2); // dw/dy=0
            }
        }
    }

    fields[0].makeSpectral();
}
}  // namespace chflow
#endif
#endif