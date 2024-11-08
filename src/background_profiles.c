#include "background_profiles.h"

int make_density_profile( gridConfiguration *gridCfg, 
                          double cntrl_para, 
                          double n_e[NX/2][Ny/2][Nz/2] ) {
//{{{
    // This function allows to defines the plasma density profiles. The
    // parameter "cntrl_para" allows to switch between varies options. The
    // plasma density is written into the array "n_e", which has half of the 
    // resoltion of the full-wave grid (half of the size into each direction),
    // as 1 Yee-cell consists of 2 grid points and we don't want the background
    // to vary inside of 1 Yee-cell.
    // Option 5 might be the most usefull option, as it simply reads the data
    // from an hdf5 file (for further details, see below in the code).

    size_t
        ii, jj, kk, 
        ne_start_z;
    double
        ne_max,
        ne_k0Ln,
        aux;

    // if density is larger than this value, FDTD code becomes instable
    ne_max  = period * 2./5.;

    if ( ne_profile == 1 ) {
        // plasma mirror
        for (ii=0 ; ii<(NX/2) ; ++ii) {
            for (jj=0 ; jj<(Ny/2) ; ++jj) {
                //for (kk=((gridCfg->Nz-gridCfg->d_absorb-4)/2) ; kk<(gridCfg->Nz/2) ; ++kk) {
                for (kk=0 ; kk<(Nz/2) ; ++kk) {
                    // z = m*y + b = delta_z/delta_y * y + z_start
                    //             = cntrl_para * y + z_start
                    if ( (double)kk > (cntrl_para*(double)jj + ((double)Nz-(double)d_absorb-4.)/2.) ) {
                        n_e[ii][jj][kk] = ne_max;
                    }
                }
            }
        }
    } else if ( ne_profile == 2 ) {
        // linearly increasing profile with k0Ln as slope
        // n_e(z) = m*z 
        // with m = 2*pi/(k0Ln*lambda) 
        //      z = z-starting_position
        ne_start_z  = (d_absorb + period)/2;
        //ne_start_z  = (gridCfg->d_absorb + 224.155)/2;    // benchmark scenario from STEP project: .15m/l_0*period
        if (ne_start_z%2 != 0)
            ne_start_z  += 1;
        ne_k0Ln     = cntrl_para;
        printf( "make_density_profile: ne_profile = %d, ne_start_z = %ld, k0Ln = %f\n", 
                ne_profile, ne_start_z, ne_k0Ln );
        for (ii=0 ; ii<(NX/2) ; ++ii) {
            for (jj=0 ; jj<(Ny/2) ; ++jj) {
                for (kk=0 ; kk<(Nz/2) ; ++kk) {
                    aux = ((double)kk - (double)ne_start_z) * (2.*M_PI / (ne_k0Ln*period));
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
    } else if ( ne_profile == 3 ) {
        // plasma cylinder (2D gauss) in center of grid
        //
        //                               (y-y0)^2        (z-z0)^2
        // n_e(y,z) = n_e,max * exp( -( ------------ + ------------- ) )
        //                               2*sigma_y^2    2*sigma_z^2
        //
        for (ii=0 ; ii<(NX/2) ; ++ii) {
            for (jj=0 ; jj<(Ny/2) ; ++jj) {
                for (kk=0 ; kk<(Nz/2) ; ++kk) {
                    n_e[ii][jj][kk]    = exp( -1.* (
                                                 pow((double)jj-(double)Ny/4., 2)/(2*pow(period/2.,2)) 
                                                +pow((double)kk-(double)Nz/4., 2)/(2*pow(period/2.,2))
                                             )) * 5.;
                }
            }
        }
    } else if ( ne_profile == 4 ) {
        // same as ne_profile = 3, but plasma cylinder is now along y
        for (ii=0 ; ii<(NX/2) ; ++ii) {
            for (jj=0 ; jj<(Ny/2) ; ++jj) {
                for (kk=0 ; kk<(Nz/2) ; ++kk) {
                    n_e[ii][jj][kk]    = exp( -1.* (
                                                 pow((double)ii-(double)NX/4., 2)/(2*pow(period/2.,2)) 
                                                +pow((double)kk-(double)Nz/4., 2)/(2*pow(period/2.,2))
                                             )) * 2;//5.;
                }
            }
        }
    } else if ( ne_profile == 5 ) {
        // filename and dataset-name are currently hard-coded
        //   ==> change this
        //       either provide additional parameter in function call
        //       or not load the profile here, but directly in main
        //readMyHDF( NX/2, Ny/2, Nz/2, "input/grid.h5", "n_e", n_e );
    }
    return EXIT_SUCCESS;
}//}}}


int make_B0_profile( gridConfiguration *gridCfg, 
                     double cntrl_para, 
                     double J_B0[NX][Ny][Nz] ) {
//{{{
    size_t
        ii, jj, kk; 
//    double
//        aux;

    // B0x: even-odd-odd
    // B0y: odd-even-odd
    // B0z: odd-odd-even
    if ( B0_profile == 1 ) {
        // constant field
        for (ii=0 ; ii<NX ; ii+=2) {
            for (jj=0 ; jj<Ny ; jj+=2) {
                for (kk=0 ; kk<Nz ; kk+=2) {
                    J_B0[ii  ][jj+1][kk+1] = cntrl_para;
                    J_B0[ii+1][jj  ][kk+1] = cntrl_para*.0;
                    J_B0[ii+1][jj+1][kk  ] = cntrl_para*.0;
                }
            }
        }
    } else if ( B0_profile == 2 ) {
        // filename and dataset-name are currently hard-coded
        // (just an urgent, quick, and dirty hack, sorry)
        //   ==> change this
        //       either provide additional parameter in function call
        //       or not load the profile here, but directly in main
        
        double (*B0_tmp)[Ny/2][Nz/2]  = calloc(NX/2, sizeof *B0_tmp);

        // B0_x
        readMyHDF( NX/2, Ny/2, Nz/2, "input/grid.h5", "B0_x", B0_tmp );
        for (ii=0 ; ii<NX ; ii+=2) {
            for (jj=0 ; jj<Ny ; jj+=2) {
                for (kk=0 ; kk<Nz ; kk+=2) {
                    J_B0[ii  ][jj+1][kk+1] = B0_tmp[ii/2][jj/2][kk/2];
                    B0_tmp[ii/2][jj/2][kk/2] = .0;
                }
            }
        }
        // B0_y
        readMyHDF( NX/2, Ny/2, Nz/2, "input/grid.h5", "B0_y", B0_tmp );
        for (ii=0 ; ii<NX ; ii+=2) {
            for (jj=0 ; jj<Ny ; jj+=2) {
                for (kk=0 ; kk<Nz ; kk+=2) {
                    J_B0[ii+1][jj  ][kk+1] = B0_tmp[ii/2][jj/2][kk/2];
                    B0_tmp[ii/2][jj/2][kk/2] = .0;
                }
            }
        }
        // B0_z
        readMyHDF( NX/2, Ny/2, Nz/2, "input/grid.h5", "B0_z", B0_tmp );
        for (ii=0 ; ii<NX ; ii+=2) {
            for (jj=0 ; jj<Ny ; jj+=2) {
                for (kk=0 ; kk<Nz ; kk+=2) {
                    J_B0[ii+1][jj+1][kk  ] = B0_tmp[ii/2][jj/2][kk/2];
                    B0_tmp[ii/2][jj/2][kk/2] = .0;
                }
            }
        }

        free( B0_tmp );
    }

    return EXIT_SUCCESS;
}//}}}
