#include "boundaries.h"
#include "focal.h"

int set_densityInAbsorber_v2( gridConfiguration *gridCfg, 
                              char absorber[], 
                              double n_e[gridCfg->Nx/2][gridCfg->Ny/2][gridCfg->Nz/2] ) {
//{{{

    double
        x, y, z,
        x0,x1, y0,y1, z0,z1,
        ne_absorb,              // electron density in absorber
        smooth,                 // 1: relatively steep, .2: more smooth
        ne_dist,                // proportional to distance to absorber, where n_e starts to decrease
        scale_fact;

    ne_absorb   = .0;
    smooth      = .5;//.2;
    ne_dist     = round( gridCfg->period/1 );

    x0          = (double)gridCfg->d_absorb + ne_dist;
    x1          = (double)gridCfg->Nx - (gridCfg->d_absorb + ne_dist);
    y0          = (double)gridCfg->d_absorb + ne_dist;
    y1          = (double)gridCfg->Ny - (gridCfg->d_absorb + ne_dist);
    z0          = (double)gridCfg->d_absorb + ne_dist;
    z1          = (double)gridCfg->Nz - (gridCfg->d_absorb + ne_dist);

    // scale to density grid which is only half the size of FDTD-wavefields grid
    // since 2 numerical grid points correspond to one "physical" grid point
    // and we want not to vary the background parameters within one physical
    // grid point
    x0  = round(x0*.5);
    x1  = round(x1*.5);
    y0  = round(y0*.5);
    y1  = round(y1*.5);
    z0  = round(z0*.5);
    z1  = round(z1*.5);

    // the string "absorber" is used to set in which absorber n_e will be modified
    // the comparison is done with the strstr() function, which return the address
    // of the substring if found, NULL otherwise

    // set density in x0 absorber
    if ( strstr(absorber,"x1") ) {
        for ( x=0; x<(gridCfg->Nx/2) ; ++x ) {
            scale_fact  = +.5*(    tanh(smooth*(x-x0)) + 1);        // x0 boundary
            //printf( "x1: x=%.1f, scale_fact=%f\n", x, scale_fact) ;
            for ( y=0. ; y<(gridCfg->Ny/2) ; ++y )   {
                for (z=0 ; z<(gridCfg->Nz/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }
    // set density in x1 absorber
    if ( strstr(absorber,"x2") ) {
        for ( x=0; x<(gridCfg->Nx/2) ; ++x ) {
            scale_fact  = +.5*(-1.*tanh(smooth*(x-x1)) + 1);       // x1 boundary
            //printf( "x2: x=%.1f, scale_fact=%f\n", x, scale_fact) ;
            for ( y=0. ; y<(gridCfg->Ny/2) ; ++y )   {
                for (z=0 ; z<(gridCfg->Nz/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }

    // set density in y0 absorber
    if ( strstr(absorber,"y1") ) {
        for ( y=0; y<(gridCfg->Ny/2) ; ++y ) {
            scale_fact  = +.5*(    tanh(smooth*(y-y0)) + 1);        // y0 boundary
            //printf( "y1: y=%.1f, scale_fact=%f\n", y, scale_fact) ;
            for ( x=0; x<(gridCfg->Nx/2) ; ++x ) {
                for (z=0 ; z<(gridCfg->Nz/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }
    // set density in y1 absorber
    if ( strstr(absorber,"y2") ) {
        for ( y=0; y<(gridCfg->Ny/2) ; ++y ) {
            scale_fact  = +.5*(-1.*tanh(smooth*(y-y1)) + 1);       // y1 boundary
            //printf( "y2: y=%.1f, scale_fact=%f\n", y, scale_fact) ;
            for ( x=0; x<(gridCfg->Nx/2) ; ++x ) {
                for (z=0 ; z<(gridCfg->Nz/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }

    // set density in z0 absorber
    if ( strstr(absorber,"z1") ) {
        for ( z=0 ; z<(gridCfg->Nz/2) ; ++z) {
            scale_fact  = +.5*(    tanh(smooth*(z-z0)) + 1);        // z0 boundary
            //printf( "z1: z=%.1f, scale_fact=%f\n", z, scale_fact) ;
            for ( x=0; x<(gridCfg->Nx/2) ; ++x ) {
                for ( y=0; y<(gridCfg->Ny/2) ; ++y ) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }
    // set density in z1 absorber
    if ( strstr(absorber,"z2") ) {
        for ( z=0 ; z<(gridCfg->Nz/2) ; ++z) {
            scale_fact  = +.5*(-1.*tanh(smooth*(z-z1)) + 1);       // z1 boundary
            //printf( "z2: z=%.1f, scale_fact=%f\n", z, scale_fact) ;
            for ( x=0; x<(gridCfg->Nx/2) ; ++x ) {
                for ( y=0; y<(gridCfg->Ny/2) ; ++y ) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }

    return EXIT_SUCCESS;
} //}}}


int apply_absorber( gridConfiguration *gridCfg, 
                    double eco, 
                    double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] ) {
//{{{
    size_t
        ii, jj, kk;
    double
        damp;

    // z1 absorber: z=0...d_absorb
//#pragma omp parallel for collapse(2) default(shared) private(k,j,damp) <-- can collapse be used here?
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->d_absorb-2 ; kk+=2) {
                damp = ((double)kk-(double)gridCfg->d_absorb)/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
//                if ((ii%10 == 0) && (jj%10 == 0) && (kk%10 == 0)) 
//                    printf( "z1: ii=%3d, jj=%3d, kk=%3d, (kk-d_abs)/d_abs=%f, damp=%f\n", 
//                            ii, jj, kk, ((double)kk-(double)d_absorb)/(double)d_absorb, damp );
            }
        }
    }
    // z2 absorber: z=d_absorb...Nz
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=(gridCfg->Nz-gridCfg->d_absorb) ; kk<gridCfg->Nz-2 ; kk+=2) {      //Nz-d_absorb-2 ???
                damp = ((double)kk-((double)gridCfg->Nz-(double)gridCfg->d_absorb))/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }      
    // x1 absorber: x=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
            for (ii=2 ; ii<gridCfg->d_absorb-2 ; ii+=2) {
                damp = ((double)ii-(double)gridCfg->d_absorb)/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // x2 absorber: x=d_absorb...Nx
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {  
            for (ii=(gridCfg->Nx-gridCfg->d_absorb) ; ii<gridCfg->Nx-2 ; ii+=2) {    //Nx-d_absorb-2 ???
                damp = ((double)ii-((double)gridCfg->Nx-(double)gridCfg->d_absorb))/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // y1 absorber: y=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
            for (jj=2 ; jj<gridCfg->d_absorb-2 ; jj+=2) {
                damp = ((double)jj-(double)gridCfg->d_absorb)/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // y2 absorber: y=d_absorb...Ny
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
            for (jj=(gridCfg->Ny-gridCfg->d_absorb) ; jj<gridCfg->Ny-2 ; jj+=2) {  //Ny-d_absorb-2 ???
                damp = ((double)jj-((double)gridCfg->Ny-(double)gridCfg->d_absorb))/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int apply_absorber_ref( gridConfiguration *gridCfg, 
                        double eco, 
                        double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref] ) {
//{{{
    size_t
        ii, jj, kk;
    double
        damp;

    // z1 absorber: z=0...d_absorb
//#pragma omp parallel for collapse(2) default(shared) private(k,j,damp) <-- can collapse be used here?
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->d_absorb-2 ; kk+=2) {
                damp = ((double)kk-(double)gridCfg->d_absorb)/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
//                if ((ii%10 == 0) && (jj%10 == 0) && (kk%10 == 0)) 
//                    printf( "z1: ii=%3d, jj=%3d, kk=%3d, (kk-d_abs)/d_abs=%f, damp=%f\n", 
//                            ii, jj, kk, ((double)kk-(double)d_absorb)/(double)d_absorb, damp );
            }
        }
    }
    // z2 absorber: z=d_absorb...Nz
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=(gridCfg->Nz_ref-gridCfg->d_absorb) ; kk<gridCfg->Nz_ref-2 ; kk+=2) {      //Nz-d_absorb-2 ???
                damp = ((double)kk-((double)gridCfg->Nz_ref-(double)gridCfg->d_absorb))/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }      
    // x1 absorber: x=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            for (ii=2 ; ii<gridCfg->d_absorb-2 ; ii+=2) {
                damp = ((double)ii-(double)gridCfg->d_absorb)/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // x2 absorber: x=d_absorb...Nx
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {  
            for (ii=(gridCfg->Nx-gridCfg->d_absorb) ; ii<gridCfg->Nx-2 ; ii+=2) {    //Nx-d_absorb-2 ???
                damp = ((double)ii-((double)gridCfg->Nx-(double)gridCfg->d_absorb))/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // y1 absorber: y=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            for (jj=2 ; jj<gridCfg->d_absorb-2 ; jj+=2) {
                damp = ((double)jj-(double)gridCfg->d_absorb)/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // y2 absorber: y=d_absorb...Ny
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            for (jj=(gridCfg->Ny-gridCfg->d_absorb) ; jj<gridCfg->Ny-2 ; jj+=2) {  //Ny-d_absorb-2 ???
                damp = ((double)jj-((double)gridCfg->Ny-(double)gridCfg->d_absorb))/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int apply_absorber_v2( size_t N_x, size_t N_y, size_t N_z, int d_absorb, double eco, 
                       char absorber[],
                       double EB_WAVE[N_x][N_y][N_z] ) {
//{{{
    size_t
        ii, jj, kk;
    double
        damp;

    // z1 absorber: z=0...d_absorb
//#pragma omp parallel for collapse(2) default(shared) private(k,j,damp) <-- can collapse be used here?
    if ( strstr(absorber,"z1") ) {      
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (jj=2 ; jj<N_y-2 ; jj+=2) {
                for (kk=2 ; kk<d_absorb-2 ; kk+=2) {
                    damp = ((double)kk-(double)d_absorb)/(double)d_absorb;
                    damp = ABSORBER_DAMPING(eco,damp);

                    EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                    EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                    EB_WAVE[ii  ][jj  ][kk+1] *= damp;
//                if ((ii%10 == 0) && (jj%10 == 0) && (kk%10 == 0)) 
//                    printf( "z1: ii=%3d, jj=%3d, kk=%3d, (kk-d_abs)/d_abs=%f, damp=%f\n", 
//                            ii, jj, kk, ((double)kk-(double)d_absorb)/(double)d_absorb, damp );
                }
            }
        }
    }
    // z2 absorber: z=d_absorb...Nz
    if ( strstr(absorber,"z2") ) {      
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (jj=2 ; jj<N_y-2 ; jj+=2) {
                for (kk=(N_z-d_absorb) ; kk<N_z-2 ; kk+=2) {      //Nz-d_absorb-2 ???
                    damp = ((double)kk-((double)N_z-(double)d_absorb))/(double)d_absorb;
                    damp = ABSORBER_DAMPING(eco,damp);

                    EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                    EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                    EB_WAVE[ii  ][jj  ][kk+1] *= damp;
                }
            }
        }      
    }
    // x1 absorber: x=0...d_absorb
    if ( strstr(absorber,"x1") ) {
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (jj=2 ; jj<N_y-2 ; jj+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {
                for (ii=2 ; ii<d_absorb-2 ; ii+=2) {
                    damp = ((double)ii-(double)d_absorb)/(double)d_absorb;
                    damp = ABSORBER_DAMPING(eco,damp);

                    EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                    EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                    EB_WAVE[ii  ][jj  ][kk+1] *= damp;
                }
            }
        }
    }
    if ( strstr(absorber,"x2") ) {
    // x2 absorber: x=d_absorb...Nx
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (jj=2 ; jj<N_y-2 ; jj+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {  
                for (ii=(N_x-d_absorb) ; ii<N_x-2 ; ii+=2) {    //Nx-d_absorb-2 ???
                    damp = ((double)ii-((double)N_x-(double)d_absorb))/(double)d_absorb;
                    damp = ABSORBER_DAMPING(eco,damp);

                    EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                    EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                    EB_WAVE[ii  ][jj  ][kk+1] *= damp;
                }
            }
        }
    }
    // y1 absorber: y=0...d_absorb
    if ( strstr(absorber,"y1") ) {
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {
                for (jj=2 ; jj<d_absorb-2 ; jj+=2) {
                    damp = ((double)jj-(double)d_absorb)/(double)d_absorb;
                    damp = ABSORBER_DAMPING(eco,damp);

                    EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                    EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                    EB_WAVE[ii  ][jj  ][kk+1] *= damp;
                }
            }
        }
    }
    // y2 absorber: y=d_absorb...Ny
    if ( strstr(absorber,"y2") ) {
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {
                for (jj=(N_y-d_absorb) ; jj<N_y-2 ; jj+=2) {  //Ny-d_absorb-2 ???
                    damp = ((double)jj-((double)N_y-(double)d_absorb))/(double)d_absorb;
                    damp = ABSORBER_DAMPING(eco,damp);

                    EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                    EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                    EB_WAVE[ii  ][jj  ][kk+1] *= damp;
                }
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int apply_numerical_viscosity( gridConfiguration *gridCfg,
                               double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] ) {
    //{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk;

    double
        aux, ny;

    ny  = 1e-4;
    
#pragma omp parallel for default(shared) private(ii,jj,kk,aux)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
                aux = ny*(   EB_WAVE[ii  ][jj  ][kk+1+2] 
                          +  EB_WAVE[ii  ][jj  ][kk+1-2]
                          -2*EB_WAVE[ii  ][jj  ][kk+1  ]
                        );
                EB_WAVE[ii  ][jj  ][kk+1] += aux;
            }
        }
    }

    return EXIT_SUCCESS;
}//}}}


int abc_Mur_saveOldE_xdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[8][gridCfg->Ny][gridCfg->Nz] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        jj, kk, 
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
            // store values at x=0 and x=1
            // Ex: odd-even-even
            E_old[0+1][jj  ][kk  ]  = EB_WAVE[0+offset+1][jj  ][kk  ];
            E_old[2+1][jj  ][kk  ]  = EB_WAVE[2+offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_old[0  ][jj+1][kk  ]  = EB_WAVE[0+offset  ][jj+1][kk  ];
            E_old[2  ][jj+1][kk  ]  = EB_WAVE[2+offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_old[0  ][jj  ][kk+1]  = EB_WAVE[0+offset  ][jj  ][kk+1];
            E_old[2  ][jj  ][kk+1]  = EB_WAVE[2+offset  ][jj  ][kk+1];

            // store values at x=Nx-1 and x=Nx-2
            // Ex: odd-even-even
            E_old[4+1][jj  ][kk  ]  = EB_WAVE[gridCfg->Nx-4-offset+1][jj  ][kk  ];
            E_old[6+1][jj  ][kk  ]  = EB_WAVE[gridCfg->Nx-2-offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_old[4  ][jj+1][kk  ]  = EB_WAVE[gridCfg->Nx-4-offset  ][jj+1][kk  ];
            E_old[6  ][jj+1][kk  ]  = EB_WAVE[gridCfg->Nx-2-offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_old[4  ][jj  ][kk+1]  = EB_WAVE[gridCfg->Nx-4-offset  ][jj  ][kk+1];
            E_old[6  ][jj  ][kk+1]  = EB_WAVE[gridCfg->Nx-2-offset  ][jj  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldE_ydir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[gridCfg->Nx][8][gridCfg->Nz] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, kk,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
            // store values at y=0 and y=1
            // Ex: odd-even-even
            E_old[ii+1][0  ][kk  ]  = EB_WAVE[ii+1][0+offset  ][kk  ];
            E_old[ii+1][2  ][kk  ]  = EB_WAVE[ii+1][2+offset  ][kk  ];
            // Ey: even-odd-even
            E_old[ii  ][0+1][kk  ]  = EB_WAVE[ii  ][0+offset+1][kk  ];
            E_old[ii  ][2+1][kk  ]  = EB_WAVE[ii  ][2+offset+1][kk  ];
            // Ez: even-even-odd
            E_old[ii  ][0  ][kk+1]  = EB_WAVE[ii  ][0+offset  ][kk+1];
            E_old[ii  ][2  ][kk+1]  = EB_WAVE[ii  ][2+offset  ][kk+1];

            // store values at x=Nx-1 and x=Nx-2
            // Ex: odd-even-even
            E_old[ii+1][4  ][kk  ]  = EB_WAVE[ii+1][gridCfg->Ny-4-offset  ][kk  ];
            E_old[ii+1][6  ][kk  ]  = EB_WAVE[ii+1][gridCfg->Ny-2-offset  ][kk  ];
            // Ey: even-odd-even
            E_old[ii  ][4+1][kk  ]  = EB_WAVE[ii  ][gridCfg->Ny-4-offset+1][kk  ];
            E_old[ii  ][6+1][kk  ]  = EB_WAVE[ii  ][gridCfg->Ny-2-offset+1][kk  ];
            // Ez: even-even-odd
            E_old[ii  ][4  ][kk+1]  = EB_WAVE[ii  ][gridCfg->Ny-4-offset  ][kk+1];
            E_old[ii  ][6  ][kk+1]  = EB_WAVE[ii  ][gridCfg->Ny-2-offset  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldE_zdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[gridCfg->Nx][gridCfg->Ny][8] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            // store values at z=0 and z=1
            // Ex: odd-even-even
            E_old[ii+1][jj  ][0  ]  = EB_WAVE[ii+1][jj  ][0+offset  ];
            E_old[ii+1][jj  ][2  ]  = EB_WAVE[ii+1][jj  ][2+offset  ];
            // Ey: even-odd-even
            E_old[ii  ][jj+1][0  ]  = EB_WAVE[ii  ][jj+1][0+offset  ];
            E_old[ii  ][jj+1][2  ]  = EB_WAVE[ii  ][jj+1][2+offset  ];
            // Ez: even-even-odd
            E_old[ii  ][jj  ][0+1]  = EB_WAVE[ii  ][jj  ][0+offset+1];
            E_old[ii  ][jj  ][2+1]  = EB_WAVE[ii  ][jj  ][2+offset+1];

            // store values at z=Nz-1 and z=Nz-2
            // Ex: odd-even-even
            E_old[ii+1][jj  ][4  ]  = EB_WAVE[ii+1][jj  ][gridCfg->Nz-4-offset  ];
            E_old[ii+1][jj  ][6  ]  = EB_WAVE[ii+1][jj  ][gridCfg->Nz-2-offset  ];
            // Ey: even-odd-even
            E_old[ii  ][jj+1][4  ]  = EB_WAVE[ii  ][jj+1][gridCfg->Nz-4-offset  ];
            E_old[ii  ][jj+1][6  ]  = EB_WAVE[ii  ][jj+1][gridCfg->Nz-2-offset  ];
            // Ez: even-even-odd
            E_old[ii  ][jj  ][4+1]  = EB_WAVE[ii  ][jj  ][gridCfg->Nz-4-offset+1];
            E_old[ii  ][jj  ][6+1]  = EB_WAVE[ii  ][jj  ][gridCfg->Nz-2-offset+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldEref_xdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[8][gridCfg->Ny][gridCfg->Nz_ref] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        jj, kk, 
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            // store values at x=0 and x=1
            // Ex: odd-even-even
            E_old[0+1][jj  ][kk  ]  = EB_WAVE_ref[0+offset+1][jj  ][kk  ];
            E_old[2+1][jj  ][kk  ]  = EB_WAVE_ref[2+offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_old[0  ][jj+1][kk  ]  = EB_WAVE_ref[0+offset  ][jj+1][kk  ];
            E_old[2  ][jj+1][kk  ]  = EB_WAVE_ref[2+offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_old[0  ][jj  ][kk+1]  = EB_WAVE_ref[0+offset  ][jj  ][kk+1];
            E_old[2  ][jj  ][kk+1]  = EB_WAVE_ref[2+offset  ][jj  ][kk+1];

            // store values at x=Nx-1 and x=Nx-2
            // Ex: odd-even-even
            E_old[4+1][jj  ][kk  ]  = EB_WAVE_ref[gridCfg->Nx-4-offset+1][jj  ][kk  ];
            E_old[6+1][jj  ][kk  ]  = EB_WAVE_ref[gridCfg->Nx-2-offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_old[4  ][jj+1][kk  ]  = EB_WAVE_ref[gridCfg->Nx-4-offset  ][jj+1][kk  ];
            E_old[6  ][jj+1][kk  ]  = EB_WAVE_ref[gridCfg->Nx-2-offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_old[4  ][jj  ][kk+1]  = EB_WAVE_ref[gridCfg->Nx-4-offset  ][jj  ][kk+1];
            E_old[6  ][jj  ][kk+1]  = EB_WAVE_ref[gridCfg->Nx-2-offset  ][jj  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldEref_ydir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[gridCfg->Nx][8][gridCfg->Nz_ref] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, kk,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            // store values at y=0 and y=1
            // Ex: odd-even-even
            E_old[ii+1][0  ][kk  ]  = EB_WAVE_ref[ii+1][0+offset  ][kk  ];
            E_old[ii+1][2  ][kk  ]  = EB_WAVE_ref[ii+1][2+offset  ][kk  ];
            // Ey: even-odd-even
            E_old[ii  ][0+1][kk  ]  = EB_WAVE_ref[ii  ][0+offset+1][kk  ];
            E_old[ii  ][2+1][kk  ]  = EB_WAVE_ref[ii  ][2+offset+1][kk  ];
            // Ez: even-even-odd
            E_old[ii  ][0  ][kk+1]  = EB_WAVE_ref[ii  ][0+offset  ][kk+1];
            E_old[ii  ][2  ][kk+1]  = EB_WAVE_ref[ii  ][2+offset  ][kk+1];

            // store values at x=Nx-1 and x=Nx-2
            // Ex: odd-even-even
            E_old[ii+1][4  ][kk  ]  = EB_WAVE_ref[ii+1][gridCfg->Ny-4-offset  ][kk  ];
            E_old[ii+1][6  ][kk  ]  = EB_WAVE_ref[ii+1][gridCfg->Ny-2-offset  ][kk  ];
            // Ey: even-odd-even
            E_old[ii  ][4+1][kk  ]  = EB_WAVE_ref[ii  ][gridCfg->Ny-4-offset+1][kk  ];
            E_old[ii  ][6+1][kk  ]  = EB_WAVE_ref[ii  ][gridCfg->Ny-2-offset+1][kk  ];
            // Ez: even-even-odd
            E_old[ii  ][4  ][kk+1]  = EB_WAVE_ref[ii  ][gridCfg->Ny-4-offset  ][kk+1];
            E_old[ii  ][6  ][kk+1]  = EB_WAVE_ref[ii  ][gridCfg->Ny-2-offset  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldEref_zdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[gridCfg->Nx][gridCfg->Ny][8] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            // store values at z=0 and z=1
            // Ex: odd-even-even
            E_old[ii+1][jj  ][0  ]  = EB_WAVE_ref[ii+1][jj  ][0+offset  ];
            E_old[ii+1][jj  ][2  ]  = EB_WAVE_ref[ii+1][jj  ][2+offset  ];
            // Ey: even-odd-even
            E_old[ii  ][jj+1][0  ]  = EB_WAVE_ref[ii  ][jj+1][0+offset  ];
            E_old[ii  ][jj+1][2  ]  = EB_WAVE_ref[ii  ][jj+1][2+offset  ];
            // Ez: even-even-odd
            E_old[ii  ][jj  ][0+1]  = EB_WAVE_ref[ii  ][jj  ][0+offset+1];
            E_old[ii  ][jj  ][2+1]  = EB_WAVE_ref[ii  ][jj  ][2+offset+1];

            // store values at z=Nz-1 and z=Nz-2
            // Ex: odd-even-even
            E_old[ii+1][jj  ][4  ]  = EB_WAVE_ref[ii+1][jj  ][gridCfg->Nz_ref-4-offset  ];
            E_old[ii+1][jj  ][6  ]  = EB_WAVE_ref[ii+1][jj  ][gridCfg->Nz_ref-2-offset  ];
            // Ey: even-odd-even
            E_old[ii  ][jj+1][4  ]  = EB_WAVE_ref[ii  ][jj+1][gridCfg->Nz_ref-4-offset  ];
            E_old[ii  ][jj+1][6  ]  = EB_WAVE_ref[ii  ][jj+1][gridCfg->Nz_ref-2-offset  ];
            // Ez: even-even-odd
            E_old[ii  ][jj  ][4+1]  = EB_WAVE_ref[ii  ][jj  ][gridCfg->Nz_ref-4-offset+1];
            E_old[ii  ][jj  ][6+1]  = EB_WAVE_ref[ii  ][jj  ][gridCfg->Nz_ref-2-offset+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_1st( gridConfiguration *gridCfg, 
                 char absorber[],
                 double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                 double E_old_xdir[8][gridCfg->Ny][gridCfg->Nz], 
                 double E_old_ydir[gridCfg->Nx][8][gridCfg->Nz], 
                 double E_old_zdir[gridCfg->Nx][gridCfg->Ny][8] ) {
//{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk,
        offset;             // refers to EB_WAVE only

    double
        cnst;

    cnst    = (gridCfg->dt-gridCfg->dx)/(gridCfg->dt+gridCfg->dx);
    offset  = 2;

    // the string "absorber" is used to set which absorber is treated
    // the comparison is done with the strstr() function, which return the address
    // of the substring if found, NULL otherwise
    // NOTE: "if (strstr(absorber,"x1))" should be sufficient
    //       "if (strstr(absorber,"x1) != NULL)" should be equivalent

    // absorber into x-direction
    if ( strstr(absorber,"x1") ) {      
        //printf("abs_Mur_1st_v2: x1\n");
#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
                // absorber at x=0 grid boundary
                // Ex: odd-even-even
                EB_WAVE[offset+0+1][jj  ][kk  ] = E_old_xdir[2+1][jj  ][kk  ]
                    + cnst * (    EB_WAVE[offset+2+1][jj  ][kk  ] 
                              -E_old_xdir[0+1       ][jj  ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[offset+0  ][jj+1][kk  ] = E_old_xdir[2  ][jj+1][kk  ]
                    + cnst * (    EB_WAVE[offset+2  ][jj+1][kk  ] 
                              -E_old_xdir[0         ][jj+1][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[offset+0  ][jj  ][kk+1] = E_old_xdir[2  ][jj  ][kk+1]
                    + cnst * (    EB_WAVE[offset+2  ][jj  ][kk+1] 
                              -E_old_xdir[0         ][jj  ][kk+1] );
            }
        }
    }
    if ( strstr(absorber,"x2") ) {
        //printf("abs_Mur_1st_v2: x2\n");
#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
                // absorber at x=Nx grid boundary
                // Ex: odd-even-even
                EB_WAVE[gridCfg->Nx-2-offset+1][jj  ][kk  ]    = E_old_xdir[4+1][jj  ][kk  ]
                    + cnst * (    EB_WAVE[gridCfg->Nx-4-offset+1][jj  ][kk  ] 
                              -E_old_xdir[6+1                   ][jj  ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[gridCfg->Nx-2-offset  ][jj+1][kk  ]    = E_old_xdir[4  ][jj+1][kk  ]
                    + cnst * (    EB_WAVE[gridCfg->Nx-4-offset  ][jj+1][kk  ] 
                              -E_old_xdir[6                     ][jj+1][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[gridCfg->Nx-2-offset  ][jj  ][kk+1]    = E_old_xdir[4  ][jj  ][kk+1]
                    + cnst * (    EB_WAVE[gridCfg->Nx-4-offset  ][jj  ][kk+1] 
                              -E_old_xdir[6                     ][jj  ][kk+1] );
            }
        }
    }

    // absorber into y-direction
    if ( strstr(absorber,"y1") ) {
        //printf("abs_Mur_1st_v2: y1\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
        for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
            for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
                // absorber at y=0 grid boundary
                // Ex: odd-even-even
                EB_WAVE[ii+1][offset+0  ][kk  ] = E_old_ydir[ii+1][2  ][kk  ]
                    + cnst * (    EB_WAVE[ii+1][offset+2  ][kk  ]
                              -E_old_ydir[ii+1][0         ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[ii  ][offset+0+1][kk  ] = E_old_ydir[ii  ][2+1][kk  ]
                    + cnst * (    EB_WAVE[ii  ][offset+2+1][kk  ]
                              -E_old_ydir[ii  ][0+1       ][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[ii  ][offset+0  ][kk+1] = E_old_ydir[ii  ][2  ][kk+1]
                    + cnst * (    EB_WAVE[ii  ][offset+2  ][kk+1]
                              -E_old_ydir[ii  ][0         ][kk+1] );
            }
        }
    }
    if ( strstr(absorber,"y2") ) {
        //printf("abs_Mur_1st_v2: y2\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
        for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
            for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
                // absorber at y=Ny grid boundary
                // Ex: odd-even-even
                EB_WAVE[ii+1][gridCfg->Ny-2-offset  ][kk  ] = E_old_ydir[ii+1][4  ][kk  ]
                    + cnst * (    EB_WAVE[ii+1][gridCfg->Ny-4-offset  ][kk  ]
                              -E_old_ydir[ii+1][6                     ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[ii  ][gridCfg->Ny-2-offset+1][kk  ] = E_old_ydir[ii  ][4+1][kk  ]
                    + cnst * (    EB_WAVE[ii  ][gridCfg->Ny-4-offset+1][kk  ]
                              -E_old_ydir[ii  ][6+1                   ][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[ii  ][gridCfg->Ny-2-offset  ][kk+1] = E_old_ydir[ii  ][4  ][kk+1]
                    + cnst * (    EB_WAVE[ii  ][gridCfg->Ny-4-offset  ][kk+1]
                              -E_old_ydir[ii  ][6                     ][kk+1] );
            }
        }
    }

    // absorber into z-direction
    if ( strstr(absorber,"z1") ) {
        //printf("abs_Mur_1st_v2: z1\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
        for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
            for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
                // absorber at z=0 grid boundary
                // Ex: odd-even-even
                EB_WAVE[ii+1][jj  ][offset+0]   = E_old_zdir[ii+1][jj  ][2  ]
                    + cnst * (    EB_WAVE[ii+1][jj  ][offset+2  ]
                              -E_old_zdir[ii+1][jj  ][0  ]        );
                // Ey: even-odd-even
                EB_WAVE[ii  ][jj+1][offset+0]   = E_old_zdir[ii  ][jj+1][2  ]
                    + cnst * (    EB_WAVE[ii  ][jj+1][offset+2  ]
                              -E_old_zdir[ii  ][jj+1][0  ]        );
                // Ez: even-even-odd
                EB_WAVE[ii  ][jj  ][offset+0+1] = E_old_zdir[ii  ][jj  ][2+1]
                    + cnst * (    EB_WAVE[ii  ][jj  ][offset+2+1]
                              -E_old_zdir[ii  ][jj  ][0+1]        );
            }
        }
    }
    if ( strstr(absorber,"z2") ) {
        //printf("abs_Mur_1st_v2: z2\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
        for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
            for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
                // absorber at z=Nz grid boundary
                // Ex: odd-even-even
                EB_WAVE[ii+1][jj  ][gridCfg->Nz-2-offset  ]    = E_old_zdir[ii+1][jj  ][4  ]
                    + cnst * (    EB_WAVE[ii+1][jj  ][gridCfg->Nz-4-offset  ]
                              -E_old_zdir[ii+1][jj  ][6  ]     );
                // Ey: even-odd-even
                EB_WAVE[ii  ][jj+1][gridCfg->Nz-2-offset  ]    = E_old_zdir[ii  ][jj+1][4  ]
                    + cnst * (    EB_WAVE[ii  ][jj+1][gridCfg->Nz-4-offset  ]
                              -E_old_zdir[ii  ][jj+1][6  ]     );
                // Ez: even-even-odd
                EB_WAVE[ii  ][jj  ][gridCfg->Nz-2-offset+1]    = E_old_zdir[ii  ][jj  ][4+1]
                    + cnst * (    EB_WAVE[ii  ][jj  ][gridCfg->Nz-4-offset+1]
                              -E_old_zdir[ii  ][jj  ][6+1]     );
            }
        }
    }

    return EXIT_SUCCESS;

} //}}}


int abc_Mur_1st_ref( gridConfiguration *gridCfg,
                     double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                     double E_old_xdir[8][gridCfg->Ny][gridCfg->Nz_ref], 
                     double E_old_ydir[gridCfg->Nx][8][gridCfg->Nz_ref], 
                     double E_old_zdir[gridCfg->Nx][gridCfg->Ny][8] ) {
//{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk,
        offset;             // refers to EB_WAVE only

    double
        cnst;

    cnst    = (gridCfg->dt-gridCfg->dx)/(gridCfg->dt+gridCfg->dx);
    offset  = 2;

    // absorber into x-direction
#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            // absorber at x=0 grid boundary
            // Ex: odd-even-even
            EB_WAVE[offset+0+1][jj  ][kk  ] = E_old_xdir[2+1][jj  ][kk  ]
                + cnst * (    EB_WAVE[offset+2+1][jj  ][kk  ] 
                          -E_old_xdir[0+1       ][jj  ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[offset+0  ][jj+1][kk  ] = E_old_xdir[2  ][jj+1][kk  ]
                + cnst * (    EB_WAVE[offset+2  ][jj+1][kk  ] 
                          -E_old_xdir[0         ][jj+1][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[offset+0  ][jj  ][kk+1] = E_old_xdir[2  ][jj  ][kk+1]
                + cnst * (    EB_WAVE[offset+2  ][jj  ][kk+1] 
                          -E_old_xdir[0         ][jj  ][kk+1] );
            // absorber at x=Nx grid boundary
            // Ex: odd-even-even
            EB_WAVE[gridCfg->Nx-2-offset+1][jj  ][kk  ]    = E_old_xdir[4+1][jj  ][kk  ]
                + cnst * (    EB_WAVE[gridCfg->Nx-4-offset+1][jj  ][kk  ] 
                          -E_old_xdir[6+1                   ][jj  ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[gridCfg->Nx-2-offset  ][jj+1][kk  ]    = E_old_xdir[4  ][jj+1][kk  ]
                + cnst * (    EB_WAVE[gridCfg->Nx-4-offset  ][jj+1][kk  ] 
                          -E_old_xdir[6                     ][jj+1][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[gridCfg->Nx-2-offset  ][jj  ][kk+1]    = E_old_xdir[4  ][jj  ][kk+1]
                + cnst * (    EB_WAVE[gridCfg->Nx-4-offset  ][jj  ][kk+1] 
                          -E_old_xdir[6                     ][jj  ][kk+1] );
        }
    }

    // absorber into y-direction
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            // absorber at y=0 grid boundary
            // Ex: odd-even-even
            EB_WAVE[ii+1][offset+0  ][kk  ] = E_old_ydir[ii+1][2  ][kk  ]
                + cnst * (    EB_WAVE[ii+1][offset+2  ][kk  ]
                          -E_old_ydir[ii+1][0         ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[ii  ][offset+0+1][kk  ] = E_old_ydir[ii  ][2+1][kk  ]
                + cnst * (    EB_WAVE[ii  ][offset+2+1][kk  ]
                          -E_old_ydir[ii  ][0+1       ][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[ii  ][offset+0  ][kk+1] = E_old_ydir[ii  ][2  ][kk+1]
                + cnst * (    EB_WAVE[ii  ][offset+2  ][kk+1]
                          -E_old_ydir[ii  ][0         ][kk+1] );
            // absorber at y=Ny grid boundary
            // Ex: odd-even-even
            EB_WAVE[ii+1][gridCfg->Ny-2-offset  ][kk  ] = E_old_ydir[ii+1][4  ][kk  ]
                + cnst * (    EB_WAVE[ii+1][gridCfg->Ny-4-offset  ][kk  ]
                          -E_old_ydir[ii+1][6                     ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[ii  ][gridCfg->Ny-2-offset+1][kk  ] = E_old_ydir[ii  ][4+1][kk  ]
                + cnst * (    EB_WAVE[ii  ][gridCfg->Ny-4-offset+1][kk  ]
                          -E_old_ydir[ii  ][6+1                   ][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[ii  ][gridCfg->Ny-2-offset  ][kk+1] = E_old_ydir[ii  ][4  ][kk+1]
                + cnst * (    EB_WAVE[ii  ][gridCfg->Ny-4-offset  ][kk+1]
                          -E_old_ydir[ii  ][6                     ][kk+1] );
        }
    }

    // absorber into z-direction
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            // absorber at z=0 grid boundary
            // Ex: odd-even-even
            EB_WAVE[ii+1][jj  ][offset+0]   = E_old_zdir[ii+1][jj  ][2  ]
                + cnst * (    EB_WAVE[ii+1][jj  ][offset+2  ]
                          -E_old_zdir[ii+1][jj  ][0  ]        );
            // Ey: even-odd-even
            EB_WAVE[ii  ][jj+1][offset+0]   = E_old_zdir[ii  ][jj+1][2  ]
                + cnst * (    EB_WAVE[ii  ][jj+1][offset+2  ]
                          -E_old_zdir[ii  ][jj+1][0  ]        );
            // Ez: even-even-odd
            EB_WAVE[ii  ][jj  ][offset+0+1] = E_old_zdir[ii  ][jj  ][2+1]
                + cnst * (    EB_WAVE[ii  ][jj  ][offset+2+1]
                          -E_old_zdir[ii  ][jj  ][0+1]        );
            // absorber at z=Nz grid boundary
            // Ex: odd-even-even
            EB_WAVE[ii+1][jj  ][gridCfg->Nz_ref-2-offset  ]    = E_old_zdir[ii+1][jj  ][4  ]
                + cnst * (    EB_WAVE[ii+1][jj  ][gridCfg->Nz_ref-4-offset  ]
                          -E_old_zdir[ii+1][jj  ][6  ]     );
            // Ey: even-odd-even
            EB_WAVE[ii  ][jj+1][gridCfg->Nz_ref-2-offset  ]    = E_old_zdir[ii  ][jj+1][4  ]
                + cnst * (    EB_WAVE[ii  ][jj+1][gridCfg->Nz_ref-4-offset  ]
                          -E_old_zdir[ii  ][jj+1][6  ]     );
            // Ez: even-even-odd
            EB_WAVE[ii  ][jj  ][gridCfg->Nz_ref-2-offset+1]    = E_old_zdir[ii  ][jj  ][4+1]
                + cnst * (    EB_WAVE[ii  ][jj  ][gridCfg->Nz_ref-4-offset+1]
                          -E_old_zdir[ii  ][jj  ][6+1]     );
        }
    }

    return EXIT_SUCCESS;

} //}}}


/*split PML related functions*/
/*void initialize_split_PML(pmlBoundary *PML, gridConfiguration *gridCfg, int pml_size){

    printf("Initializing Split-fields for PML Boundary conditions. \n");

    ALLOC_3D(PML->EBx,  gridCfg->Nx, gridCfg->Ny, gridCfg->Nz, double);
    ALLOC_3D(PML->EBy,  gridCfg->Nx, gridCfg->Ny, gridCfg->Nz, double);
    ALLOC_3D(PML->EBz,  gridCfg->Nx, gridCfg->Ny, gridCfg->Nz, double);
    ALLOC_3D(PML->EB,   gridCfg->Nx, gridCfg->Ny, gridCfg->Nz, double);

    ALLOC_1D(PML->bx,  gridCfg->Nx, double);
    ALLOC_1D(PML->by,  gridCfg->Ny, double);
    ALLOC_1D(PML->bz,  gridCfg->Nz, double);
    ALLOC_1D(PML->cx,  gridCfg->Nx, double);
    ALLOC_1D(PML->cy,  gridCfg->Ny, double);
    ALLOC_1D(PML->cz,  gridCfg->Nz, double);

    split_PML_parameter(PML, gridCfg, pml_size);
    
}

void apply_split_PML( pmlBoundary *PML, gridConfiguration *gridCfg, int pml_size,
                      double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] ){

    int ii, jj, kk;

    #pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for ( ii=2 ; ii<gridCfg->Nx-2 ; ii+=2 ) {
        for ( jj=2 ; jj<gridCfg->Ny-2 ; jj+=2 ) {
            for ( kk=2 ; kk<gridCfg->Nz-2 ; kk+=2 ) {
                if(ii <= pml_size + 2 || ii >= gridCfg->Nx - 2 - pml_size ||
                   jj <= pml_size + 2 || jj >= gridCfg->Ny - 2 - pml_size || 
                   kk <= pml_size + 2 || kk >= gridCfg->Nz - 2 - pml_size ){

                    if(ii == pml_size + 2 && jj == pml_size + 2 && kk == pml_size + 2){
                        printf("EBxy = %f \n", EBx(ii ,jj+1 ,kk) );
                        printf("EBxz = %f \n", EBx(ii ,jj ,kk+1) );
                        printf("EB_wave = %f \n", EB_WAVE[ii+1][jj+1][kk] );
                        printf("EB_wave = %f \n", EB_WAVE[ii+1][jj-1][kk] );
                    }

                    //Ex split fields
                    EBx(ii ,jj+1 ,kk) =  ( by(jj+1)*EBx(ii ,jj+1 ,kk) + 
                                           cz(kk+1)*(EB_WAVE[ii+1][jj+1][kk] - EB_WAVE[ii+1][jj-1][kk]) );

                    EBx(ii ,jj ,kk+1) =  ( bz(kk+1)*EBx(ii ,jj ,kk+1) - 
                                           cy(jj+1)*(EB_WAVE[ii+1][jj][kk+1] - EB_WAVE[ii+1][jj][kk-1]) );

                    //Ey split fields
                    EBy(ii+1 ,jj ,kk) =  ( bx(ii+1)*EBy(ii+1 ,jj ,kk) - 
                                           cz(kk+1)*(EB_WAVE[ii+1][jj+1][kk] - EB_WAVE[ii-1][jj+1][kk]) );

                    EBy(ii ,jj ,kk+1) =  ( bz(kk+1)*EBy(ii ,jj ,kk+1) + 
                                           cx(ii+1)*(EB_WAVE[ii][jj+1][kk+1] - EB_WAVE[ii][jj+1][kk-1]) );

                    //Ez split fields
                    EBz(ii+1 ,jj ,kk) =  ( bx(ii+1)*EBz(ii+1 ,jj ,kk) + 
                                           cy(jj+1)*(EB_WAVE[ii+1][jj][kk+1] - EB_WAVE[ii-1][jj][kk+1]) );

                    EBz(ii ,jj+1 ,kk) =  ( by(jj+1)*EBz(ii ,jj+1 ,kk) - 
                                           cx(ii+1)*(EB_WAVE[ii][jj+1][kk+1] - EB_WAVE[ii][jj-1][kk+1]) );

                    //Bx split fields
                    EBx(ii+1 ,jj ,kk+1) =  ( by(jj)*EBx(ii+1 ,jj ,kk+1) + 
                                             cz(kk)*(EB_WAVE[ii][jj+2][kk+1] - EB_WAVE[ii][jj][kk+1]) );

                    EBx(ii+1 ,jj+1 ,kk) =  ( bz(kk)*EBx(ii+1 ,jj+1 ,kk) - 
                                             cy(jj)*(EB_WAVE[ii][jj+1][kk+2] - EB_WAVE[ii][jj+1][kk]) );

                    //By split fields
                    EBy(ii ,jj+1 ,kk+1) =  ( bx(ii)*EBy(ii ,jj+1 ,kk+1) - 
                                             cz(kk)*(EB_WAVE[ii+2][jj][kk+1] - EB_WAVE[ii][jj][kk+1]) );

                    EBy(ii+1 ,jj+1 ,kk) =  ( bz(kk)*EBy(ii+1 ,jj+1 ,kk) + 
                                             cx(ii)*(EB_WAVE[ii+1][jj][kk+2] - EB_WAVE[ii+1][jj][kk]) );

                    //Bz split fields
                    EBz(ii ,jj+1 ,kk+1) =  ( bx(ii)*EBz(ii ,jj+1 ,kk+1) + 
                                             cy(jj)*(EB_WAVE[ii+2][jj+1][kk] - EB_WAVE[ii][jj+1][kk]) );

                    EBz(ii+1 ,jj ,kk+1) =  ( by(jj)*EBz(ii+1 ,jj ,kk+1) - 
                                             cx(ii)*(EB_WAVE[ii+1][jj+2][kk] - EB_WAVE[ii+1][jj][kk]) );
                }
            }
        }
    }
}

void update_split_PML(  pmlBoundary *PML, gridConfiguration *gridCfg, int pml_size, 
                        double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz]){
    
    int ii, jj, kk;

    #pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for ( ii=2 ; ii<gridCfg->Nx-2 ; ii+=2 ) {
        for ( jj=2 ; jj<gridCfg->Ny-2 ; jj+=2 ) {
            for ( kk=2 ; kk<gridCfg->Nz-2 ; kk+=2 ) {
                if(ii <= pml_size + 2 || ii >= gridCfg->Nx - 2 - pml_size ||
                   jj <= pml_size + 2 || jj >= gridCfg->Ny - 2 - pml_size || 
                   kk <= pml_size + 2 || kk >= gridCfg->Nz - 2 - pml_size){

                    //Electric field update
                    EB_WAVE[ii+1][jj][kk] = EBx(ii ,jj+1 ,kk) + EBx(ii ,jj ,kk+1);
                    EB_WAVE[ii][jj+1][kk] = EBy(ii+1 ,jj ,kk) + EBy(ii ,jj ,kk+1);
                    EB_WAVE[ii][jj][kk+1] = EBz(ii+1 ,jj ,kk) + EBz(ii ,jj+1 ,kk);

                    //Magnetic field update
                    EB_WAVE[ii][jj+1][kk+1] = EBx(ii+1 ,jj ,kk+1) + EBx(ii+1 ,jj+1 ,kk);
                    EB_WAVE[ii+1][jj][kk+1] = EBy(ii ,jj+1 ,kk+1) + EBy(ii+1 ,jj+1 ,kk);
                    EB_WAVE[ii+1][jj+1][kk] = EBz(ii ,jj+1 ,kk+1) + EBz(ii+1 ,jj ,kk+1);
                   }
            }
        }
    }
}

double sigma(int pml_size, double nn, int m, double ds){

    double sig, sig_max, R_0;

    R_0 = pow(10,-6);
    sig_max = -(m+1)*log( R_0 )/(2*pml_size*ds);

    sig = pow( (nn) /(pml_size), m) * sig_max;

    return sig;  
}


void split_PML_parameter(pmlBoundary *PML, gridConfiguration *gridCfg, int pml_size){
    
    int ii, jj, kk, count;
    double sig, c;
    
    count = pml_size;
    c = gridCfg->dt/gridCfg->dx;

    //#pragma omp parallel  
    for ( ii=2 ; ii<gridCfg->Nx-2 ; ii+=2 ) {
        if(ii <= pml_size + 2){

            sig = sigma(pml_size, count - 1, 4, gridCfg->dx);
            //Electric component
            bx(ii+1) = (2 - sig*gridCfg->dt)/(2 + sig*gridCfg->dt);
            cx(ii+1) = c * ( 2/(2 + sig*gridCfg->dt) );
            
            sig = sigma(pml_size, count, 4, gridCfg->dx);
            //Magnetic component
            bx(ii) = (2 - sig*gridCfg->dt)/(2 + sig*gridCfg->dt);
            cx(ii) = c * ( 2/(2 + sig*gridCfg->dt) );
            
            //printf("Hcx = %f \n", cx(ii) );
            //printf("Ecx = %f \n", cx(ii+1) );

            count -= 2;
        }else if( ii >= gridCfg->Nx - 2 - pml_size){
            count += 2; 

            sig = sigma(pml_size, count - 1, 4, gridCfg->dx);
            //Electric component
            bx(ii+1) = (2 - sig*gridCfg->dt)/(2 + sig*gridCfg->dt);
            cx(ii+1) = c * ( 2/(2 + sig*gridCfg->dt) );
            
            sig = sigma(pml_size, count, 4, gridCfg->dx);
            //Magnetic component
            bx(ii) = (2 - sig*gridCfg->dt)/(2 + sig*gridCfg->dt);
            cx(ii) = c * ( 2/(2 + sig*gridCfg->dt) );  

            //printf("Hcx = %f \n", cx(ii) );
            //printf("Ecx = %f \n", cx(ii+1) );  
        }
    }

    for ( jj=2 ; jj<gridCfg->Ny-2 ; jj+=2 ) {
        if(jj <= pml_size + 2){

            sig = sigma(pml_size, count - 1, 4, gridCfg->dx);
            //Electric component
            by(jj+1) = (2 - sig*gridCfg->dt)/(2 + sig*gridCfg->dt);
            cy(jj+1) = c * ( 2/(2 + sig*gridCfg->dt) );

            sig = sigma(pml_size, count, 4, gridCfg->dx);
            //Magnetic component
            by(jj) = (2 - sig*gridCfg->dt)/(2 + sig*gridCfg->dt);
            cy(jj) = c * ( 2/(2 + sig*gridCfg->dt) );    

            count -= 2;
        }else if( jj >= gridCfg->Ny - 2 - pml_size){
            count += 2;

            sig = sigma(pml_size, count - 1, 4, gridCfg->dx);
            //Electric component
            by(jj+1) = (2 - sig*gridCfg->dt)/(2 + sig*gridCfg->dt);
            cy(jj+1) = c * ( 2/(2 + sig*gridCfg->dt) );
            
            sig = sigma(pml_size, count, 4, gridCfg->dx);
            //Magnetic component
            by(jj) = (2 - sig*gridCfg->dt)/(2 + sig*gridCfg->dt);
            cy(jj) = c * ( 2/(2 + sig*gridCfg->dt) );    
        }
    }
    
    for ( kk=2 ; kk<gridCfg->Nz-2 ; kk+=2 ) {
        if(kk <= pml_size + 2){

            sig = sigma(pml_size, count - 1, 4, gridCfg->dx);
            //Electric component
            bz(kk+1) = (2 - sig*gridCfg->dt)/(2 + sig*gridCfg->dt);
            cz(kk+1) = c * ( 2/(2 + sig*gridCfg->dt) );
            
            sig = sigma(pml_size, count, 4, gridCfg->dx);
            //Magnetic component
            bz(kk) = (2 - sig*gridCfg->dt)/(2 + sig*gridCfg->dt);
            cz(kk) = c * ( 2/(2 + sig*gridCfg->dt) );    

            count -= 2;

        }else if( kk >= gridCfg->Nz - 2 - pml_size){

            count += 2; 

            sig = sigma(pml_size, count - 1, 4, gridCfg->dx);
            //Electric component
            bz(kk+1) = (2 - sig*gridCfg->dt)/(2 + sig*gridCfg->dt);
            cz(kk+1) = c * ( 2/(2 + sig*gridCfg->dt) );
            
            sig = sigma(pml_size, count, 4, gridCfg->dx);
            //Magnetic component
            bz(kk) = (2 - sig*gridCfg->dt)/(2 + sig*gridCfg->dt);
            cz(kk) = c * ( 2/(2 + sig*gridCfg->dt) );    
        }
    }
}*/


/*Apply boundary conditions to simulation*/
void setBoundary(   gridConfiguration *gridCfg,   
                    namePath *pathFile,
                    abcBoundary *ABC){

    int pml_size;
    pml_size = gridCfg->d_absorb;

    if(pathFile->boundary == 1){ //normal absorbing boundary condition
        printf("Absorbing Boundary conditions applied. \n");
        ABC->eco = 10./(double)(gridCfg->period);
    }

    if(pathFile->boundary == 2){
        //initialize_split_PML(PML, gridCfg);
    }    

    /*if(pathFile->boundary == 3){
        initialize_split_PML(PML, gridCfg, pml_size);
    }*/

}

void computeBoundary(   gridConfiguration *gridCfg,   
                        namePath *pathFile, 
                        abcBoundary *ABC,
                        double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], int t_int,
                        double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref]){

    int pml_size;
    pml_size = gridCfg->d_absorb;

    if(pathFile->boundary == 1){ //normal absorbing boundary condition
        
        /*printf("No problem here \n");
        apply_absorber(     &gridCfg, eco, EB_WAVE );
        printf("No problem here \n");
        apply_absorber_ref( &gridCfg, eco, EB_WAVE_ref );*/

    }

    if(pathFile->boundary == 2){
        
    }    

    /*if(pathFile->boundary == 3){
        //initialize_split_PML(PML, gridCfg);
        apply_split_PML( PML, gridCfg, pml_size, EB_WAVE);
        update_split_PML( PML, gridCfg, pml_size, EB_WAVE);
        
    }*/

}