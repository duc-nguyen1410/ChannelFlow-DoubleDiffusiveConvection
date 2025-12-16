/**
 * 
 * Original author: Duc Nguyen
 */

 #include <iomanip>
 #include <iostream>
 
 
 #include <sys/stat.h>
 #include <fstream>
 #include <iomanip>
 #include <iostream>
 #include <string>
 
 #include "cfbasics/cfvector.h"
 #include "cfbasics/mathdefs.h"
 #include "channelflow/chebyshev.h"
 #include "channelflow/dns.h"
 #include "channelflow/flowfield.h"
 #include "channelflow/poissonsolver.h"
 #include "channelflow/symmetry.h"
 #include "channelflow/tausolver.h"
 #include "channelflow/utilfuncs.h"
 #include "modules/ddc/macros.h"
 #include "modules/ddc/ddcdsi.h"
 #include "modules/ddc/ddc.h"
 #include "modules/ddc/boundaryCondition.h"
 #include "modules/ddc/turbulenceStatistics.h"
 using namespace std;
 using namespace chflow;
 
 string printdiagnostics(FlowField& u, const DNS& dns, Real t, const TimeStep& dt, Real nu, Real umin, bool vardt,
                         bool pl2norm, bool pchnorm, bool pdissip, bool pshear, bool pdiverge, bool pUbulk, bool pubulk,
                         bool pdPdx, bool pcfl);
 
 int main(int argc, char* argv[]) {
     cfMPI_Init(&argc, &argv);
     {
         WriteProcessInfo(argc, argv);
         string purpose(
             "integrate wall-bounded double-diffusive convection "
             "from a given initial condition and save fields to disk.");
 
         ArgList args(argc, argv, purpose);
 
         DDCFlags flags(args);
         TimeStep dt(flags);
 
         args.section("Program options");
 
         const string outdir = args.getpath("-o", "--outdir", "data/", "output directory");
         const string ulabel = args.getstr("-ul", "--ulabel", "u", "output velocity field prefix");
         const string tlabel = args.getstr("-tl", "--tlabel", "t", "output temperature field prefix");
         const string slabel = args.getstr("-sl", "--slabel", "s", "output salinity field prefix");
 
         const bool pcfl = args.getflag("-cfl", "--cfl", "print CFL number each dT");
         const bool pl2norm = args.getflag("-l2", "--l2norm", "print L2Norm(u) each dT");
         const bool pchnorm = args.getbool("-ch", "--chebyNorm", true, "print chebyNorm(u) each dT");
         const bool pdissip = args.getflag("-D", "--dissipation", "print dissipation each dT");
         const bool pshear = args.getflag("-I", "--input", "print wall shear power input each dT");
         const bool pdiverge = args.getflag("-dv", "--divergence", "print divergence each dT");
         const bool pubulk = args.getflag("-u", "--ubulk", "print ubulk each dT");
         const bool pUbulk = args.getflag("-Up", "--Ubulk-print", "print Ubulk each dT");
         const bool pdPdx = args.getflag("-p", "--pressure", "print pressure gradient each dT");
         const Real umin = args.getreal("-u", "--umin", 0.0, "stop if chebyNorm(u) < umin");
 
         const Real ecfmin = args.getreal("-e", "--ecfmin", 0.0, "stop if Ecf(u) < ecfmin");
         const int cadence = args.getint("-c", "--cadence", 10, "recompute CFL per c time steps");
         const int saveint = args.getint("-s", "--saveinterval", 1, "save fields every s dT");
         
         const int nproc0 =
             args.getint("-np0", "--nproc0", 0, "number of MPI-processes for transpose/number of parallel ffts");
         const int nproc1 = args.getint("-np1", "--nproc1", 0, "number of MPI-processes for one fft");
         
 
         const string outdir_profiles = args.getpath("-op", "--outdir_profiles", "profiles/", "output directory of mean profiles");
         const bool savetot = args.getflag("-savetot", "--savetotfields", "save total fields instead of perturbations");
         const bool average = args.getflag("-average", "--averagefields", "average fields");
         
         const bool saveslicex = args.getflag("-saveslicex", "--saveSliceX", "save profile at slice x");
         const int slicexidx = args.getint("-slicex", "--sliceXIndex", 0,"index of slice x_i, 0<=i<Nx");
         const bool saveslicey = args.getflag("-saveslicey", "--saveSliceY", "save profile at slice y");
         const int sliceyidx = args.getint("-slicey", "--sliceYIndex", 0,"index of slice y_j, 0<=j<Ny");
         
         const string uname = args.getstr(3, "<flowfield>", "initial guess for the velocity solution");
         const string tname = args.getstr(2, "<flowfield>", "initial guess for the temperature solution");
         const string sname = args.getstr(1, "<flowfield>", "initial guess for the salinity solution");
          
         args.check();
         args.save("./");
         mkdir(outdir);
         args.save(outdir);
         flags.save(outdir);
         
         // use input fields
         
         CfMPI* cfmpi = &CfMPI::getInstance(nproc0, nproc1);
 
         printout("Constructing u,q, and optimizing FFTW...");
         FlowField u(uname, cfmpi);
         FlowField temp(tname, cfmpi);
         FlowField salt(sname, cfmpi);
         
         const int Nx = u.Nx();
         const int Ny = u.Ny();
         const int Nz = u.Nz();
         const Real Lx = u.Lx();
         const Real Lz = u.Lz();
         const Real a = u.a();
         const Real b = u.b();
 
         FlowField q(Nx, Ny, Nz, 1, Lx, Lz, a, b, cfmpi); // for pressure
         const bool inttime =
             (abs(saveint * dt.dT() - int(saveint * dt.dT())) < 1e-12) && (abs(flags.t0 - int(flags.t0)) < 1e-12)
                 ? true
                 : false;
         
 
         cout << "Parameters: " << endl;
         cout << "Re = " << flags.Rey << endl;
         cout << "Ra = " << flags.Ra << endl;
         cout << "Le = " << flags.Le << endl;
         cout << "Pr = " << flags.Pr << endl;
         cout << "Ri = " << flags.Ri << endl;
         cout << "Rrho = " << flags.Rrho << endl;
         #ifdef P7
         cout << "Rsep = " << flags.Rsep << endl;
         #endif
         cout << "Ua = " << flags.ulowerwall << endl;
         cout << "Wa = " << flags.wlowerwall << endl;
         cout << "Ub = " << flags.uupperwall << endl;
         cout << "Wb = " << flags.wupperwall << endl;
         #ifdef P5
         cout << "Ta = " << flags.tlowerwall << endl;
         cout << "Tb = " << flags.tupperwall << endl;
         #endif
         #ifdef P6
         cout << "Sa = " << flags.slowerwall << endl;
         cout << "Sb = " << flags.supperwall << endl;
         #endif
 
         // Construct data fields: 3d velocity, 1d temperature, 1d salinity and 1d pressure
         vector<FlowField> fields = {u, temp, salt, q};
 
         // time-averaged fields
         FlowField aveu;
         FlowField avetemp;
         FlowField avesalt;
         int avecount = 0;
         if (average){
             aveu = u;
             avetemp = temp;
             avesalt = salt;
         }
         // vector<FlowField> avefields = {aveu, avetemp, avesalt};// time-avereged fields
 
         // Construct Navier-Stoke integrator, set integration method
         cout << "Building DDC-DNS..." << flush;
         DDC ddc(fields, flags);
         cout << "done" << endl;
 
         ddc.Ubase().save(outdir + "Ubase");
         ddc.Wbase().save(outdir + "Wbase");
         ddc.Tbase().save(outdir + "Tbase");
         ddc.Sbase().save(outdir + "Sbase");
 
         PressureSolver psolver(u, ddc.Ubase(), ddc.Wbase(), flags.nu, flags.Vsuck,
                                flags.nonlinearity);  // NOT CORRECT FOR DDC
         psolver.solve(q, u);
 
         ios::openmode openflag = (flags.t0 > 0) ? ios::app : ios::out;
         ofstream eout, x0out;
         openfile(eout, outdir + "energy.asc", openflag);
         eout << ddcfieldstatsheader_t("t", flags) << endl;
         
         #ifdef P6
         #ifdef SAVESTATS
         DDCTurbStats meanProfiles(Nx,Ny);
         mkdir(outdir_profiles);
         #endif
         #endif
 
         #ifdef FREEZEvelocity
         FlowField u0 = u;
         #endif
         int count=floorf(flags.t0/dt.dT());
         for (Real t = flags.t0; t <= flags.T; t += cadence*dt.dt()) {
             // check dt
             if(floorf(t/dt.dT())>=count){
                 
                 string s;
                 // s = printdiagnostics(fields[0], ddc, t, dt, flags.nu, umin, dt.variable(), pl2norm, pchnorm, pdissip,
                 //                     pshear, pdiverge, pUbulk, pubulk, pdPdx, pcfl);
                 // if (ecfmin > 0 && Ecf(fields[0]) < ecfmin) {
                 //     cferror("Ecf < ecfmin == " + r2s(ecfmin) + ", exiting");
                 // }
                 // cout << s;
 
                 s = ddcfieldstats_t(fields[0], fields[1], fields[2], t, flags);
                 eout << s << endl;
 
                 #ifdef P6
                 #ifdef SAVESTATS
                 meanProfiles.addSnapshot(fields, flags, saveslicex, slicexidx, saveslicey, sliceyidx);
                 meanProfiles.saveTurbStats(outdir_profiles+"meanprofile"+i2s(int(count)), outdir_profiles+"sliceY"+i2s(int(count)), fields, saveslicex, slicexidx, saveslicey, sliceyidx);
                 #endif
                 #endif
 
                 if (average){
                     aveu = 1.0/(avecount+1)*(avecount*aveu+fields[0]);
                     avetemp = 1.0/(avecount+1)*(avecount*avetemp+fields[1]);
                     avesalt = 1.0/(avecount+1)*(avecount*avesalt+fields[2]);
                     avecount+=1;
                 }
                 
                 // Write fields to disk
                 if (saveint != 0 && count % saveint == 0) {
                     if(!savetot){
                         fields[0].save(outdir + ulabel + t2s(t, inttime));//<<--- save perturbations
                         #ifdef P5
                         fields[1].save(outdir + tlabel + t2s(t, inttime));
                         #endif
                         #ifdef P6
                         fields[2].save(outdir + slabel + t2s(t, inttime));
                         #endif
 
                         // save time-averaged fields
                         if (average){
                             aveu.save(outdir + "ave"+ ulabel);//<<--- save perturbations
                             avetemp.save(outdir + "ave"+ tlabel);//<<--- save perturbations
                             avesalt.save(outdir + "ave"+ slabel);//<<--- save perturbations
                         }
                     }else{
                         FlowField u_tot = totalVelocity(fields[0], flags); u_tot.save(outdir + ulabel + t2s(t, inttime));//<<--- save total fields
                         #ifdef P5
                         FlowField temp_tot = totalTemperature(fields[1], flags); temp_tot.save(outdir + tlabel + t2s(t, inttime));
                         #endif
                         #ifdef P6
                         FlowField salt_tot = totalSalinity(fields[2], flags); salt_tot.save(outdir + slabel + t2s(t, inttime));
                         #endif
                     }
                 }
                 count+=1;
 
                 cout << "Time = "<< t << ", dt = "<< dt.dt() << ", CFL = " << ddc.CFL(fields[0]) << endl;
             }
             
             #ifdef FREEZEvelocity
             for (int step = 0; step < cadence; ++step) {
                 fields[0] = u0;
                 cout << "." << flush;
                 ddc.advance(fields, 1);
             }  
             #else
             // Take 'cadence' steps
             ddc.advance(fields, cadence);
             #endif
 
             if (dt.variable() && dt.adjust(ddc.CFL(fields[0]))) {
                 ddc.reset_dt(dt);
             }
         }
     }
     cfMPI_Finalize();
 }
 
 
 string printdiagnostics(FlowField& u, const DNS& dns, Real t, const TimeStep& dt, Real nu, Real umin, bool vardt,
                         bool pl2norm, bool pchnorm, bool pdissip, bool pshear, bool pdiverge, bool pUbulk, bool pubulk,
                         bool pdPdx, bool pcfl) {
     // Printing diagnostics
     stringstream sout;
     sout << "           t == " << t << endl;
     if (vardt)
         sout << "          dt == " << Real(dt) << endl;
     if (pl2norm)
         sout << "   L2Norm(u) == " << L2Norm(u) << endl;
 
     if (pchnorm || umin != 0.0) {
         Real chnorm = chebyNorm(u);
         sout << "chebyNorm(u) == " << chnorm << endl;
         if (chnorm < umin) {
             cout << "Exiting: chebyNorm(u) < umin." << endl;
             exit(0);
         }
     }
     Real h = 0.5 * (u.b() - u.a());
     u += dns.Ubase();
     if (pl2norm)
         sout << "   energy(u+U) == " << 0.5 * L2Norm(u) << endl;
     if (pdissip)
         sout << "   dissip(u+U) == " << dissipation(u) << endl;
     if (pshear)
         sout << "wallshear(u+U) == " << abs(wallshearLower(u)) + abs(wallshearUpper(u)) << endl;
     if (pdiverge)
         sout << "  divNorm(u+U) == " << divNorm(u) << endl;
     if (pUbulk)
         sout << "mean u+U Ubulk == " << dns.Ubulk() << endl;
     u -= dns.Ubase();
     if (u.taskid() == u.task_coeff(0, 0)) {
         if (pubulk)
             sout << "         ubulk == " << Re(u.profile(0, 0, 0)).mean() << endl;
     }
     if (pdPdx)
         sout << "          dPdx == " << dns.dPdx() << endl;
     if (pl2norm)
         sout << "     L2Norm(u) == " << L2Norm(u) << endl;
     if (pl2norm)
         sout << "   L2Norm3d(u) == " << L2Norm3d(u) << endl;
 
     Real cfl = dns.CFL(u);
     if (u.taskid() == u.task_coeff(0, 0)) {
         ChebyCoeff U = dns.Ubase();
         ChebyCoeff W = dns.Wbase();
 
         U.makeSpectral();
         U += Re(u.profile(0, 0, 0));
         Real Ucenter = U.eval(0.5 * (u.a() + u.b()));
         Real Uwall = pythag(0.5 * (U.eval_b() - U.eval_a()), 0.5 * (W.eval_b() - W.eval_a()));
         Real Umean = U.mean();
         // sout << "        1/nu == " << 1 / nu << endl;
         sout << "  Uwall h == " << Uwall * h << endl;
         sout << "  Ubulk h == " << dns.Ubulk() * h << endl;
         sout << "  Umean h == " << Umean << " * " << h << endl;
         sout << "  Umean h == " << Umean * h << endl;
         sout << " Uparab h == " << 1.5 * dns.Ubulk() * h << endl;
         sout << "Ucenter h == " << Ucenter * h << endl;
     }
     sout << "         CFL == " << cfl << endl;
     return sout.str();
 }