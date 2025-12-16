
#include <fstream>
#include <iomanip>
#include <iostream>
#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/symmetry.h"
#include "channelflow/utilfuncs.h"
#include "cfbasics/mathdefs.h"
#include "modules/ddc/macros.h"
using namespace std;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        string purpose("Reload a field and save it");
        ArgList args(argc, argv, purpose);        
        const string inname = args.getstr(2, "<fieldname>", "input file");
        const string outname = args.getstr(1, "<fieldname>", "input file");
        
        args.check();
        args.save("./");

        // get input field
        FlowField u(inname);
        lint Nx = u.Nx();
        lint Ny = u.Ny();
        lint Nz = u.Nz();
        lint Nd = u.Nd();
        lint nxlocmin = u.nxlocmin();
        lint nxlocmax = u.nxlocmin() + u.Nxloc();
        lint nylocmin = u.nylocmin();
        lint nylocmax = u.nylocmax();
        u.makePhysical();

        

        // prepare input dataset
        // double data[Nx][Ny][Nz][Nd];
        std::vector<double> data(Nx * Ny * Nz * Nd, 0.0);

    #ifdef HAVE_MPI
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (lint nz = 0; nz < Nz; ++nz)
                for (lint ny = nylocmin; ny < nylocmax; ++ny) {
                    for (lint id = 0; id < u.Nd(); ++id){
                        // data[nx][ny][nz][id] = u(nx, ny, nz, id);
                        data[nx * Ny * Nz * Nd + ny * Nz * Nd + nz * Nd + id] = u(nx, ny, nz, id);
                    }
                }
    #else
        for (lint ny = nylocmin; ny < nylocmax; ++ny)
            for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
                for (lint nz = 0; nz < Nz; ++nz) {
                    for (lint id = 0; id < u.Nd(); ++id){
                        data[nx * Ny * Nz * Nd + ny * Nz * Nd + nz * Nd + id] = u(nx, ny, nz, id);
                    }
                }
    #endif

        // prepare output field
        FlowField uout(Nx, Ny, Nz, Nd, u.Lx(), u.Lz(), u.a(), u.b());
        Nz = uout.Nz();
        nxlocmin = uout.nxlocmin();
        nxlocmax = uout.nxlocmin() + uout.Nxloc();
        nylocmin = uout.nylocmin();
        nylocmax = uout.nylocmax();
        uout.makePhysical();

    #ifdef HAVE_MPI
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (lint nz = 0; nz < Nz; ++nz)
                for (lint ny = nylocmin; ny < nylocmax; ++ny) {
                    for (lint id = 0; id < Nd; ++id){
                        uout(nx, ny, nz, id) = data[nx * Ny * Nz * Nd + ny * Nz * Nd + nz * Nd + id];
                    }
                }
    #else
        for (lint ny = nylocmin; ny < nylocmax; ++ny)
            for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
                for (lint nz = 0; nz < Nz; ++nz) {
                    for (lint id = 0; id < Nd; ++id){
                        uout(nx, ny, nz, id) = data[nx * Ny * Nz * Nd + ny * Nz * Nd + nz * Nd + id];
                    }
                }
    #endif

        uout.save(outname);
    }
    cfMPI_Finalize();
}
