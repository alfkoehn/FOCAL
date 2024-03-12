#include <stdlib.h>

#include "focal.h"
#include "background_profiles.h"

int make_density_profile( gridConfiguration *gridCfg, 
                          double cntrl_para, 
                          double n_e[gridCfg->Nx/2][gridCfg->Ny/2][gridCfg->Nz/2] ) {
//{{{
    size_t
        ii, jj, kk, 
        ne_start_z;
    double
        ne_max,
        ne_k0Ln,
        aux;

    // if density is larger than this value, FDTD code becomes instable
    ne_max  = gridCfg->period * 2./5.;

    if ( gridCfg->ne_profile == 1 ) {
        // plasma mirror
        for (ii=0 ; ii<(gridCfg->Nx/2) ; ++ii) {
            for (jj=0 ; jj<(gridCfg->Ny/2) ; ++jj) {
                //for (kk=((gridCfg->Nz-gridCfg->d_absorb-4)/2) ; kk<(gridCfg->Nz/2) ; ++kk) {
                for (kk=0 ; kk<(gridCfg->Nz/2) ; ++kk) {
                    // z = m*y + b = delta_z/delta_y * y + z_start
                    //             = cntrl_para * y + z_start
                    if ( (double)kk > (cntrl_para*(double)jj + ((double)gridCfg->Nz-(double)gridCfg->d_absorb-4.)/2.) ) {
                        n_e[ii][jj][kk] = ne_max;
                    }
                }
            }
        }
    } else if ( gridCfg->ne_profile == 2 ) {
        // linearly increasing profile with k0Ln as slope
        // n_e(z) = m*z 
        // with m = 2*pi/(k0Ln*lambda) 
        //      z = z-starting_position
        ne_start_z  = (gridCfg->d_absorb + gridCfg->period)/2;
        //ne_start_z  = (gridCfg->d_absorb + 224.155)/2;    // benchmark scenario from STEP project: .15m/l_0*period
        if (ne_start_z%2 != 0)
            ne_start_z  += 1;
        ne_k0Ln     = cntrl_para;
        printf( "make_density_profile: ne_profile = %d, ne_start_z = %ld, k0Ln = %f\n", 
                gridCfg->ne_profile, ne_start_z, ne_k0Ln );
        for (ii=0 ; ii<(gridCfg->Nx/2) ; ++ii) {
            for (jj=0 ; jj<(gridCfg->Ny/2) ; ++jj) {
                for (kk=0 ; kk<(gridCfg->Nz/2) ; ++kk) {
                    aux = ((double)kk - (double)ne_start_z) * (2.*M_PI / (ne_k0Ln*gridCfg->period));
                    // negative density values are unphysical
                    if (aux < .0)
                        aux = .0;
                    // FDTD full-wave codes only allow for a maximum density value
                    if (aux > ne_max) {
                        aux  = ne_max;
                        printf( "    maximum density achieved (ii, jj, kk = %ld, %ld, %ld): %f\n", ii, jj, kk, aux );
                    }
                    //if ((ii%10 == 0) && (jj%10 == 0) && (kk%10 == 0))
                    //if (kk%10 == 0)
                    //    printf( "    n_e[%4ld][%4ld][%4ld] = %f\n", ii, jj, kk, aux );
                    n_e[ii][jj][kk] = aux;
                }
            }
        }
    } else if ( gridCfg->ne_profile == 3 ) {
        // plasma cylinder (2D gauss) in center of grid
        //
        //                               (y-y0)^2        (z-z0)^2
        // n_e(y,z) = n_e,max * exp( -( ------------ + ------------- ) )
        //                               2*sigma_y^2    2*sigma_z^2
        //
        for (ii=0 ; ii<(gridCfg->Nx/2) ; ++ii) {
            for (jj=0 ; jj<(gridCfg->Ny/2) ; ++jj) {
                for (kk=0 ; kk<(gridCfg->Nz/2) ; ++kk) {
                    n_e[ii][jj][kk]    = exp( -1.* (
                                                 pow((double)jj-(double)gridCfg->Ny/4., 2)/(2*pow(gridCfg->period/2.,2)) 
                                                +pow((double)kk-(double)gridCfg->Nz/4., 2)/(2*pow(gridCfg->period/2.,2))
                                             )) * 5.;
                }
            }
        }
    } else if ( gridCfg->ne_profile == 4 ) {
        // same as ne_profile = 3, but plasma cylinder is now along y
        for (ii=0 ; ii<(gridCfg->Nx/2) ; ++ii) {
            for (jj=0 ; jj<(gridCfg->Ny/2) ; ++jj) {
                for (kk=0 ; kk<(gridCfg->Nz/2) ; ++kk) {
                    n_e[ii][jj][kk]    = exp( -1.* (
                                                 pow((double)ii-(double)gridCfg->Nx/4., 2)/(2*pow(gridCfg->period/2.,2)) 
                                                +pow((double)kk-(double)gridCfg->Nz/4., 2)/(2*pow(gridCfg->period/2.,2))
                                             )) * 2;//5.;
                }
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int make_B0_profile( gridConfiguration *gridCfg, 
                     double cntrl_para, 
                     double J_B0[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] ) {
//{{{
    size_t
        ii, jj, kk; 
//    double
//        aux;

    // B0x: even-odd-odd
    // B0y: odd-even-odd
    // B0z: odd-odd-even
    if ( gridCfg->B0_profile == 1 ) {
        // constant field
        for (ii=0 ; ii<gridCfg->Nx ; ii+=2) {
            for (jj=0 ; jj<gridCfg->Ny ; jj+=2) {
                for (kk=0 ; kk<gridCfg->Nz ; kk+=2) {
                    J_B0[ii  ][jj+1][kk+1] = cntrl_para;
                    J_B0[ii+1][jj  ][kk+1] = cntrl_para*.0;
                    J_B0[ii+1][jj+1][kk  ] = cntrl_para*.0;
                }
            }
        }
    } 
    return EXIT_SUCCESS;
}//}}}

