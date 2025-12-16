/**
 * This file is similar to diffop functition in the original Channelflow but plus option of computing Nusselt number
 */
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "channelflow/chebyshev.h"
#include "channelflow/diffops.h"
#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/symmetry.h"
#include "channelflow/utilfuncs.h"
#include "modules/ddc/ddc.h"
// #include "modules/ddc/dde.h"
// #include "modules/ddc/ddcflags.h"
using namespace std;

using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        string purpose("compute Nusselt number for a given scalar's perturbation field |-ddy<totalT>h_y=a|");

        ArgList args(argc, argv, purpose);
        const string scalar = args.getstr("-scalar", "--scalarname", "temp", "choose \"temp\" or \"salt\"");
        const string scalarname = args.getstr(1, "<flowfield>", "input field");
        

        DDCFlags flags(args);
        
        FlowField scal(scalarname);
        scal.makeSpectral();
        if (scalar=="temp") {
            ChebyCoeff base = linearTemperatureProfile(scal.a(), scal.b(), scal.Ny(), flags);
            scal += base;// plus base flow

            // implicit method (this is wrong)
            // ComplexChebyCoeff h_scal = scal.profile(scal.Mx(),scal.Mz(),0);// <scalar>h
            // ComplexChebyCoeff ddy_hscalar = diff(h_scal);// ydiff(<scalar>h)
            // Real Nusselt = ddy_hscalar.a(); // ydiff(<scalar>h)|y=a
            // printf("Nusselt number |-ddy<totalT>h_y=a| = %f\n",Nusselt);

            // explicit method via dudy_a function
            Real Nu = abs(scal.dudy_a());
            if (scal.taskid()==0){
                printf("Nusselt number: |-ddy(<totalT>h)|_y=a| = %f\n",Nu);
                const string outname = scalarname + string("_Nu");
                ofstream os((outname+string(".asc")).c_str());
                os << setprecision(REAL_DIGITS) << Nu << endl;
                os.close();
            }
        }
        else if (scalar=="salt"){
            ChebyCoeff base = linearSalinityProfile(scal.a(), scal.b(), scal.Ny(), flags);
            scal += base;

            // explicit method via dudy_a function
            Real Sh = abs(scal.dudy_a());
            if (scal.taskid()==0){
                printf("Sherwood number: |-ddy(<totalS>h)|_y=a| = %f\n",Sh);
                const string outname = scalarname + string("_Sh");
                ofstream os((outname+string(".asc")).c_str());
                os << setprecision(REAL_DIGITS) << Sh << endl;
                os.close();
            }
            
            
        }
        else
            cferror("nusselt: please choose an option \"temp\" or \"salt\", rerun with -h option");

    }
    cfMPI_Finalize();
    return 0;
}
