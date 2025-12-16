/**
 * turbulence statistics for turbulent flow
 *
 * Original author: Duc Nguyen
 */

 #ifndef TUR_STATS_H
 #define TUR_STATS_H
 
 #include "channelflow/turbstats.h"
 #include "channelflow/flowfield.h"
 #include "modules/ddc/macros.h"
 #ifdef P5
 #ifdef SAVESTATS
 namespace chflow {
     class DDCTurbStats : public TurbStats {
         public:
         DDCTurbStats(const int Nx, const int Ny): 
             temp_m(Ny),
             temp_grad(Ny),
             temp_flux(Ny),
             salt_m(Ny),
             salt_grad(Ny),
             salt_flux(Ny),
             ui(Nx),// xz slice
             vi(Nx),
             ti(Nx),
             si(Nx),
             uj(Ny),// yz slice
             vj(Ny),
             tj(Ny),
             sj(Ny) { };
         void saveTurbStats(const string& filebase, const string& filebase_sliceY, vector<FlowField>& fieldsn, bool saveSliceX=false, const int sliceXidx=0, bool saveSliceY=false, const int sliceYidx=0){
             string filename(filebase);
             filename += string(".csv");
             ofstream os(filename.c_str());
             os << setprecision(REAL_DIGITS);
             char s = ',';
             int Ny = temp_m.length();
             Vector y = fieldsn[0].ygridpts();
             
             os  << "ypoints"
                 << s << "temp_m" << s << "temp_grad" << s << "temp_flux" 
                 #ifdef P6 
                 << s << "salt_m" <<  s << "salt_grad"<< s << "salt_flux" 
                 #endif
                 ;
             if(saveSliceX){
                 os  << s << "u_i" << s << "v_i" << s << "temp_i" << s << "salt_i";
             }
             
             os  << '\n' ;
                  
             for (int ny = 0; ny < Ny; ++ny) {
                 os  << y[ny] 
                     << s << temp_m[ny] << s << temp_grad[ny] << s << temp_flux[ny]
                     #ifdef P6
                     << s << salt_m[ny] << s << salt_grad[ny] << s << salt_flux[ny]
                     #endif
                     ;
                 if(saveSliceX){
                 os  << s << uj[ny] << s << vj[ny] << s << tj[ny] << s << sj[ny];
                 }
                 os  << '\n';
             }
             os.close();


            if(saveSliceY){
                int Nx = ui.length();
                Vector x = fieldsn[0].xgridpts();
                ofstream fb((filebase_sliceY+string(".csv")).c_str());
                fb << setprecision(REAL_DIGITS);

                fb  << "xpoints" << s << "ui" << s << "vi" << s << "ti" << s << "si" << "\n";

                for (int nx = 0; nx < Nx; ++nx) {
                    fb  << x[nx] 
                        << s << ui[nx] << s << vi[nx]
                        << s << ti[nx] 
                        #ifdef P6
                        << s << si[nx] 
                        #endif
                        << '\n';
                }

                fb.close();
            }
         };
         void addSnapshot(vector<FlowField>& fieldsn, const DDCFlags flags, bool saveSliceX=false, const int sliceXidx=0, bool saveSliceY=false, const int sliceYidx=0){
             FlowField u = totalVelocity(fieldsn[0], flags);
             FlowField t = totalTemperature(fieldsn[1], flags);
             FlowField dtdy = ydiff(t);
             u.makePhysical();
             t.makePhysical();
             dtdy.makePhysical();
             #ifdef P6
             FlowField s = totalSalinity(fieldsn[2], flags);
             FlowField dsdy = ydiff(s);
             s.makePhysical();
             dsdy.makePhysical();
             #endif
             lint Nx = u.Nx();
             lint Ny = u.Ny();
             lint Nz = u.Nz();
             lint nxlocmin = u.nxlocmin();
             lint nxlocmax = u.nxlocmin() + u.Nxloc();
             lint nylocmin = u.nylocmin();
             lint nylocmax = u.nylocmax();
             Vector uslice(Ny),vslice(Ny),tslice(Ny),sslice(Ny);
             uj.setToZero();vj.setToZero();tj.setToZero();sj.setToZero();
             
             for (lint ny = 0; ny < Ny; ++ny){
                 #ifdef P6
                 Vector sum(6);
                 #else
                 Vector sum(3);
                 #endif
                 if (nylocmin <= ny && ny < nylocmax){
                     for (lint nx = nxlocmin; nx < nxlocmax; ++nx){
                         for (lint nz = 0; nz < Nz; ++nz) {
                             sum[0] += t(nx,ny,nz,0);
                             sum[1] += dtdy(nx,ny,nz,0);
                             sum[2] += u(nx,ny,nz,1)*t(nx,ny,nz,0);
                             #ifdef P6
                             sum[3] += s(nx,ny,nz,0);
                             sum[4] += dsdy(nx,ny,nz,0);
                             sum[5] += u(nx,ny,nz,1)*s(nx,ny,nz,0);
                             #endif
                         }
                     }
 
                     if (saveSliceX){
                        int idx = sliceXidx;
                        if (nxlocmin <= idx && idx < nxlocmax){
                            uslice[ny] = u(idx,ny,0,0);
                            vslice[ny] = u(idx,ny,0,1);
                            tslice[ny] = t(idx,ny,0,0);
                            sslice[ny] = s(idx,ny,0,0);
                        }
                         
                    }
                }
                 // Sum up results from all processes    
                 #ifdef HAVE_MPI
                 Vector tmp(sum) ;
                 MPI_Allreduce(&tmp[0], &sum[0], 1, MPI_DOUBLE, MPI_SUM, *u.comm_world());
                 MPI_Allreduce(&tmp[1], &sum[1], 1, MPI_DOUBLE, MPI_SUM, *u.comm_world());
                 MPI_Allreduce(&tmp[2], &sum[2], 1, MPI_DOUBLE, MPI_SUM, *u.comm_world());
                 #ifdef P6
                 MPI_Allreduce(&tmp[3], &sum[3], 1, MPI_DOUBLE, MPI_SUM, *u.comm_world());
                 MPI_Allreduce(&tmp[4], &sum[4], 1, MPI_DOUBLE, MPI_SUM, *u.comm_world());
                 MPI_Allreduce(&tmp[5], &sum[5], 1, MPI_DOUBLE, MPI_SUM, *u.comm_world());
                 #endif
                 #endif
                 temp_m[ny] = sum[0]/(Nx*Nz);
                 temp_grad[ny] = sum[1]/(Nx*Nz);
                 temp_flux[ny] = sum[2]/(Nx*Nz);
                 #ifdef P6
                 salt_m[ny] = sum[3]/(Nx*Nz);
                 salt_grad[ny] = sum[4]/(Nx*Nz);
                 salt_flux[ny] = sum[5]/(Nx*Nz);
                 #endif
 
                 
                #ifdef HAVE_MPI
                if (saveSliceX){
                MPI_Allreduce(&uslice[ny],&uj[ny],1,MPI_DOUBLE,MPI_SUM,*u.comm_world());
                MPI_Allreduce(&vslice[ny],&vj[ny],1,MPI_DOUBLE,MPI_SUM,*u.comm_world());
                MPI_Allreduce(&tslice[ny],&tj[ny],1,MPI_DOUBLE,MPI_SUM,*u.comm_world());
                MPI_Allreduce(&sslice[ny],&sj[ny],1,MPI_DOUBLE,MPI_SUM,*u.comm_world());
                }
                #endif
             }

             if (saveSliceY){
                int idx = sliceYidx;
                Vector usliceY(Nx),vsliceY(Nx),tsliceY(Nx),ssliceY(Nx);
                ui.setToZero();vi.setToZero();ti.setToZero();si.setToZero();
                for (lint nx = 0; nx < Nx; ++nx){
                    if (nxlocmin <= nx && nx < nxlocmax){
                        if (nylocmin <= idx && idx < nylocmax){
                            usliceY[nx] = u(nx,idx,0,0);
                            vsliceY[nx] = u(nx,idx,0,1);
                            tsliceY[nx] = t(nx,idx,0,0);
                            ssliceY[nx] = s(nx,idx,0,0);
                        }
                    }
                    #ifdef HAVE_MPI
                    MPI_Allreduce(&usliceY[nx],&ui[nx],1,MPI_DOUBLE,MPI_SUM,*u.comm_world());
                    MPI_Allreduce(&vsliceY[nx],&vi[nx],1,MPI_DOUBLE,MPI_SUM,*u.comm_world());
                    MPI_Allreduce(&tsliceY[nx],&ti[nx],1,MPI_DOUBLE,MPI_SUM,*u.comm_world());
                    MPI_Allreduce(&ssliceY[nx],&si[nx],1,MPI_DOUBLE,MPI_SUM,*u.comm_world());
                    #endif
                }
            }
         };

        
 
         private:
         Vector temp_m;// horizontally averaged profile of temperature [<T>h]
         Vector temp_grad;// horizontally averaged profile of temperature gradient [<dT/dy>h]
         Vector temp_flux;// horizontally averaged profile of temperature's convective flux [<v*T>h]
 
         Vector salt_m;// horizontally averaged profile of salinity [<S>h]
         Vector salt_grad;// horizontally averaged profile of salinity gradient [<dS/dy>h]
         Vector salt_flux;// horizontally averaged profile of salinity's convective flux [<v*S>h]
 
         // Vector dens_m;// horizontally averaged profile of density [<((\Lambda*S-T)/(\Lambda-1))>h] with \Lambda = 1/Rrho, This code use Rrho as Lambda
         // horizontally averaged profiles of turbulent diffusivity [<v*T>h/<dT/dy>h] and [<v*S>h/<dS/dy>h]
         Vector ui,vi,ti,si;
         Vector uj,vj,tj,sj;
     };
 }  // namespace chflow
 #endif
 #endif
 #endif