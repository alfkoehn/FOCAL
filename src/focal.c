#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "focal.h"

int advance_J( gridConfiguration *gridCfg, 
               double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
               double J_B0[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz],
               double n_e[gridCfg->Nx/2][gridCfg->Ny/2][gridCfg->Nz/2] ) { 
//{{{
    size_t
        ii, jj, kk;

    // Jx: odd-even-even
    // Jy: even-odd-even
    // Jz: even-even-odd
    // B0x: even-odd-odd
    // B0y: odd-even-odd
    // B0z: odd-odd-even
//#pragma omp parallel for collapse(2) default(shared) private(k,j, Jx_tmp1,Jy_tmp1,Jz_tmp1, Jx_tmp2,Jy_tmp2,Jz_tmp2, Jx_tmp3,Jy_tmp3,Jz_tmp3 ) // collapse ???
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
                // Jx: odd-even-even
                J_B0[ii+1][jj  ][kk  ]    += + gridCfg->dt*(
                        pow(2*M_PI,2) * n_e[(ii/2)][(jj/2)][(kk/2)] * EB_WAVE[ii+1][jj  ][kk  ]
                        - 2*M_PI * ( J_B0[ii  ][jj+1][kk  ]*J_B0[ii+1][jj+1][kk  ]        // +Jy*B0z
                                    -J_B0[ii  ][jj  ][kk+1]*J_B0[ii+1][jj  ][kk+1]        // -Jz*B0y
                              )
                        );
                // Jy: even-odd-even
                J_B0[ii  ][jj+1][kk  ]    += + gridCfg->dt*(
                        pow(2*M_PI,2) * n_e[(ii/2)][(jj/2)][(kk/2)] * EB_WAVE[ii  ][jj+1][kk  ]
                        -2*M_PI * (-J_B0[ii+1][jj  ][kk  ]*J_B0[ii+1][jj+1][kk  ]         // -Jx*B0z
                                   +J_B0[ii  ][jj  ][kk+1]*J_B0[ii  ][jj+1][kk+1]         // +Jz*B0x
                              )
                        );
                // Jz: even-even-odd
                J_B0[ii  ][jj  ][kk+1]    += + gridCfg->dt*(
                        pow(2*M_PI,2) * n_e[(ii/2)][(jj/2)][(kk/2)] * EB_WAVE[ii  ][jj  ][kk+1]
                        -2*M_PI * ( J_B0[ii+1][jj  ][kk  ]*J_B0[ii+1][jj  ][kk+1]        // +Jx*B0y
                                   -J_B0[ii  ][jj+1][kk  ]*J_B0[ii  ][jj+1][kk+1]        // -Jy*B0x
                              )
                        );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_B( gridConfiguration *gridCfg, 
               double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] ) {
//{{{
    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                EB_WAVE[ii  ][jj+1][kk+1]   += -1.*gridCfg->dt/gridCfg->dx*(
                        +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                        );
                // -dBy/dt = dEx/dz - dEz/dx
                EB_WAVE[ii+1][jj  ][kk+1] += -1.*gridCfg->dt/gridCfg->dx*(
                        +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                        -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                        );
                // -dBz/dt = dEy/dx - dEx/dy
                EB_WAVE[ii+1][jj+1][kk  ] += -1.*gridCfg->dt/gridCfg->dx*(
                        +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                        -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                        );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_B_ref( gridConfiguration *gridCfg, 
                   double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref] ) {
//{{{
    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                EB_WAVE[ii  ][jj+1][kk+1]   += -1.*gridCfg->dt/gridCfg->dx*(
                        +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                        );
                // -dBy/dt = dEx/dz - dEz/dx
                EB_WAVE[ii+1][jj  ][kk+1] += -1.*gridCfg->dt/gridCfg->dx*(
                        +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                        -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                        );
                // -dBz/dt = dEy/dx - dEx/dy
                EB_WAVE[ii+1][jj+1][kk  ] += -1.*gridCfg->dt/gridCfg->dx*(
                        +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                        -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                        );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_E( gridConfiguration *gridCfg, 
               double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
               double J_B0[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] ) {
//{{{
    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                EB_WAVE[ii+1][jj  ][kk  ] += gridCfg->dt/gridCfg->dx*(
                        +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                        -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                        ) - gridCfg->dt*J_B0[ii+1][jj  ][kk  ];
                // dEy/dt = (dBx/dz - dBz/dx)
                EB_WAVE[ii  ][jj+1][kk  ] += gridCfg->dt/gridCfg->dx*(
                        +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                        -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                        ) - gridCfg->dt*J_B0[ii  ][jj+1][kk  ];
                // dEz/dt = (dBy/dx - dBx/dy)
                EB_WAVE[ii  ][jj  ][kk+1] += gridCfg->dt/gridCfg->dx*(
                        +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                        ) - gridCfg->dt*J_B0[ii  ][jj  ][kk+1];
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_E_ref( gridConfiguration *gridCfg, 
                   double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref] ) { 
//{{{
    // same as advance_E but for reference fields (directional coupler)
    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                EB_WAVE[ii+1][jj  ][kk  ] += gridCfg->dt/gridCfg->dx*(
                        +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                        -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                        );
                // dEy/dt = (dBx/dz - dBz/dx)
                EB_WAVE[ii  ][jj+1][kk  ] += gridCfg->dt/gridCfg->dx*(
                        +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                        -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                        );
                // dEz/dt = (dBy/dx - dBx/dy)
                EB_WAVE[ii  ][jj  ][kk+1] += gridCfg->dt/gridCfg->dx*(
                        +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                        );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


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


