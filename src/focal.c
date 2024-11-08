#include "focal.h"

int advance_J( gridConfiguration *gridCfg, 
               double EB_WAVE[NX][NY][NZ], 
               double J_B0[NX][NY][NZ],
               double n_e[NX/2][NY/2][NZ/2] ) { 
//{{{
    // This functions advances the current density J in time. J is calculated
    // from the fluid equation of motion of the electrons and reads
    // J_new = J_old + epsion_0*w_pe^2*E - w_ce*(Jx\hat(B)_0) - nu*J
    // Note that w_pe^2 --> n_e and w_ce --> B_0 with \hat(B) being the unit
    // vector pointing into the direction of B_0.
    // nu is a term corresponding to collisional damping that can be used to
    // respresent the effect of collisional and/or avoid some numerical 
    // instabilities that might arise at resonance like the upper-hybrid
    // resonance.

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
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (jj=2 ; jj<NY-2 ; jj+=2) {
            for (kk=2 ; kk<NZ-2 ; kk+=2) {
                // Jx: odd-even-even
                J_B0[ii+1][jj  ][kk  ]    += + DT*(
                        pow(2*M_PI,2) * n_e[(ii/2)][(jj/2)][(kk/2)] * EB_WAVE[ii+1][jj  ][kk  ]
                        - 2*M_PI * ( J_B0[ii  ][jj+1][kk  ]*J_B0[ii+1][jj+1][kk  ]        // +Jy*B0z
                                    -J_B0[ii  ][jj  ][kk+1]*J_B0[ii+1][jj  ][kk+1]        // -Jz*B0y
                              )
                        );
                // Jy: even-odd-even
                J_B0[ii  ][jj+1][kk  ]    += + DT*(
                        pow(2*M_PI,2) * n_e[(ii/2)][(jj/2)][(kk/2)] * EB_WAVE[ii  ][jj+1][kk  ]
                        -2*M_PI * (-J_B0[ii+1][jj  ][kk  ]*J_B0[ii+1][jj+1][kk  ]         // -Jx*B0z
                                   +J_B0[ii  ][jj  ][kk+1]*J_B0[ii  ][jj+1][kk+1]         // +Jz*B0x
                              )
                        );
                // Jz: even-even-odd
                J_B0[ii  ][jj  ][kk+1]    += + DT*(
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
               double EB_WAVE[NX][NY][NZ] ) {
//{{{
    // B_new = B_old - nabla x E

    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (jj=2 ; jj<NY-2 ; jj+=2) {
            for (kk=2 ; kk<NZ-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                EB_WAVE[ii  ][jj+1][kk+1]   += -1.*DT/dx*(
                        +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                        );
                // -dBy/DT = dEx/dz - dEz/dx
                EB_WAVE[ii+1][jj  ][kk+1] += -1.*DT/dx*(
                        +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                        -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                        );
                // -dBz/DT = dEy/dx - dEx/dy
                EB_WAVE[ii+1][jj+1][kk  ] += -1.*DT/dx*(
                        +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                        -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                        );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_B_ref( gridConfiguration *gridCfg, 
                   double EB_WAVE[NX][NY][NZ_REF] ) {
//{{{
    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (jj=2 ; jj<NY-2 ; jj+=2) {
            for (kk=2 ; kk<NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                EB_WAVE[ii  ][jj+1][kk+1]   += -1.*DT/dx*(
                        +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                        );
                // -dBy/DT = dEx/dz - dEz/dx
                EB_WAVE[ii+1][jj  ][kk+1] += -1.*DT/dx*(
                        +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                        -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                        );
                // -dBz/DT = dEy/dx - dEx/dy
                EB_WAVE[ii+1][jj+1][kk  ] += -1.*DT/dx*(
                        +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                        -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                        );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_E( gridConfiguration *gridCfg, 
               double EB_WAVE[NX][NY][NZ], 
               double J_B0[NX][NY][NZ] ) {
//{{{
    // E_new = E_old + c^2*nablaxB - 1/epsilon_0*J

    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (jj=2 ; jj<NY-2 ; jj+=2) {
            for (kk=2 ; kk<NZ-2 ; kk+=2) {
                // dEx/DT = (dBz/dy - dBy/dz)
                EB_WAVE[ii+1][jj  ][kk  ] += DT/dx*(
                        +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                        -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                        ) - DT*J_B0[ii+1][jj  ][kk  ];
                // dEy/DT = (dBx/dz - dBz/dx)
                EB_WAVE[ii  ][jj+1][kk  ] += DT/dx*(
                        +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                        -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                        ) - DT*J_B0[ii  ][jj+1][kk  ];
                // dEz/dt = (dBy/dx - dBx/dy)
                EB_WAVE[ii  ][jj  ][kk+1] += DT/dx*(
                        +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                        ) - DT*J_B0[ii  ][jj  ][kk+1];
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_E_ref( gridConfiguration *gridCfg, 
                   double EB_WAVE[NX][NY][NZ_REF] ) { 
//{{{
    // same as advance_E but for reference fields (directional coupler)
    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (jj=2 ; jj<NY-2 ; jj+=2) {
            for (kk=2 ; kk<NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                EB_WAVE[ii+1][jj  ][kk  ] += DT/dx*(
                        +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                        -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                        );
                // dEy/dt = (dBx/dz - dBz/dx)
                EB_WAVE[ii  ][jj+1][kk  ] += DT/dx*(
                        +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                        -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                        );
                // dEz/dt = (dBy/dx - dBx/dy)
                EB_WAVE[ii  ][jj  ][kk+1] += DT/dx*(
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
                              double n_e[NX/2][NY/2][NZ/2] ) {
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
    ne_dist     = round( period/1 );

    x0          = (double)d_absorb + ne_dist;
    x1          = (double)NX - (d_absorb + ne_dist);
    y0          = (double)d_absorb + ne_dist;
    y1          = (double)NY - (d_absorb + ne_dist);
    z0          = (double)d_absorb + ne_dist;
    z1          = (double)NZ - (d_absorb + ne_dist);

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
        for ( x=0; x<(NX/2) ; ++x ) {
            scale_fact  = +.5*(    tanh(smooth*(x-x0)) + 1);        // x0 boundary
            //printf( "x1: x=%.1f, scale_fact=%f\n", x, scale_fact) ;
            for ( y=0. ; y<(NY/2) ; ++y )   {
                for (z=0 ; z<(NZ/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }
    // set density in x1 absorber
    if ( strstr(absorber,"x2") ) {
        for ( x=0; x<(NX/2) ; ++x ) {
            scale_fact  = +.5*(-1.*tanh(smooth*(x-x1)) + 1);       // x1 boundary
            //printf( "x2: x=%.1f, scale_fact=%f\n", x, scale_fact) ;
            for ( y=0. ; y<(NY/2) ; ++y )   {
                for (z=0 ; z<(NZ/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }

    // set density in y0 absorber
    if ( strstr(absorber,"y1") ) {
        for ( y=0; y<(NY/2) ; ++y ) {
            scale_fact  = +.5*(    tanh(smooth*(y-y0)) + 1);        // y0 boundary
            //printf( "y1: y=%.1f, scale_fact=%f\n", y, scale_fact) ;
            for ( x=0; x<(NX/2) ; ++x ) {
                for (z=0 ; z<(NZ/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }
    // set density in y1 absorber
    if ( strstr(absorber,"y2") ) {
        for ( y=0; y<(NY/2) ; ++y ) {
            scale_fact  = +.5*(-1.*tanh(smooth*(y-y1)) + 1);       // y1 boundary
            //printf( "y2: y=%.1f, scale_fact=%f\n", y, scale_fact) ;
            for ( x=0; x<(NX/2) ; ++x ) {
                for (z=0 ; z<(NZ/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }

    // set density in z0 absorber
    if ( strstr(absorber,"z1") ) {
        for ( z=0 ; z<(NZ/2) ; ++z) {
            scale_fact  = +.5*(    tanh(smooth*(z-z0)) + 1);        // z0 boundary
            //printf( "z1: z=%.1f, scale_fact=%f\n", z, scale_fact) ;
            for ( x=0; x<(NX/2) ; ++x ) {
                for ( y=0; y<(NY/2) ; ++y ) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }
    // set density in z1 absorber
    if ( strstr(absorber,"z2") ) {
        for ( z=0 ; z<(NZ/2) ; ++z) {
            scale_fact  = +.5*(-1.*tanh(smooth*(z-z1)) + 1);       // z1 boundary
            //printf( "z2: z=%.1f, scale_fact=%f\n", z, scale_fact) ;
            for ( x=0; x<(NX/2) ; ++x ) {
                for ( y=0; y<(NY/2) ; ++y ) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }

    return EXIT_SUCCESS;
} //}}}


int apply_absorber( gridConfiguration *gridCfg, 
                    double eco, 
                    double EB_WAVE[NX][NY][NZ] ) {
//{{{
    size_t
        ii, jj, kk;
    double
        damp;

    // z1 absorber: z=0...d_absorb
//#pragma omp parallel for collapse(2) default(shared) private(k,j,damp) <-- can collapse be used here?
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (jj=2 ; jj<NY-2 ; jj+=2) {
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
    // z2 absorber: z=d_absorb...NZ
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (jj=2 ; jj<NY-2 ; jj+=2) {
            for (kk=(NZ - d_absorb) ; kk<NZ-2 ; kk+=2) {      //NZ-d_absorb-2 ???
                damp = ((double)kk-((double)NZ-(double)d_absorb))/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }      
    // x1 absorber: x=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj<NY-2 ; jj+=2) {
        for (kk=2 ; kk<NZ-2 ; kk+=2) {
            for (ii=2 ; ii<d_absorb-2 ; ii+=2) {
                damp = ((double)ii-(double)d_absorb)/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // x2 absorber: x=d_absorb...NX
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj<NY-2 ; jj+=2) {
        for (kk=2 ; kk<NZ-2 ; kk+=2) {  
            for (ii=(NX-d_absorb) ; ii<NX-2 ; ii+=2) {    //NX-d_absorb-2 ???
                damp = ((double)ii-((double)NX-(double)d_absorb))/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // y1 absorber: y=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (kk=2 ; kk<NZ-2 ; kk+=2) {
            for (jj=2 ; jj<d_absorb-2 ; jj+=2) {
                damp = ((double)jj-(double)d_absorb)/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // y2 absorber: y=d_absorb...NY
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (kk=2 ; kk<NZ-2 ; kk+=2) {
            for (jj=(NY - d_absorb) ; jj<NY-2 ; jj+=2) {  //NY-d_absorb-2 ???
                damp = ((double)jj-((double)NY-(double)d_absorb))/(double)d_absorb;
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
                        double EB_WAVE[NX][NY][NZ_REF] ) {
//{{{
    size_t
        ii, jj, kk;
    double
        damp;

    // z1 absorber: z=0...d_absorb
//#pragma omp parallel for collapse(2) default(shared) private(k,j,damp) <-- can collapse be used here?
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (jj=2 ; jj<NY-2 ; jj+=2) {
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
    // z2 absorber: z=d_absorb...NZ
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (jj=2 ; jj<NY-2 ; jj+=2) {
            for (kk=(NZ_REF-d_absorb) ; kk<NZ_REF-2 ; kk+=2) {      //NZ-d_absorb-2 ???
                damp = ((double)kk-((double)NZ_REF-(double)d_absorb))/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }      
    // x1 absorber: x=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj<NY-2 ; jj+=2) {
        for (kk=2 ; kk<NZ_REF-2 ; kk+=2) {
            for (ii=2 ; ii<d_absorb-2 ; ii+=2) {
                damp = ((double)ii-(double)d_absorb)/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // x2 absorber: x=d_absorb...NX
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj<NY-2 ; jj+=2) {
        for (kk=2 ; kk<NZ_REF-2 ; kk+=2) {  
            for (ii=(NX-d_absorb) ; ii<NX-2 ; ii+=2) {    //NX-d_absorb-2 ???
                damp = ((double)ii-((double)NX-(double)d_absorb))/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // y1 absorber: y=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (kk=2 ; kk<NZ_REF-2 ; kk+=2) {
            for (jj=2 ; jj<d_absorb-2 ; jj+=2) {
                damp = ((double)jj-(double)d_absorb)/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // y2 absorber: y=d_absorb...NY
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (kk=2 ; kk<NZ_REF-2 ; kk+=2) {
            for (jj=(NY-d_absorb) ; jj<NY-2 ; jj+=2) {  //NY-d_absorb-2 ???
                damp = ((double)jj-((double)NY-(double)d_absorb))/(double)d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int apply_absorber_v2( size_t N_x, size_t N_y, size_t N_z, int D_absorb, double eco, 
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
                for (kk=2 ; kk<D_absorb-2 ; kk+=2) {
                    damp = ((double)kk-(double)D_absorb)/(double)D_absorb;
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
    // z2 absorber: z=d_absorb...NZ
    if ( strstr(absorber,"z2") ) {      
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (jj=2 ; jj<N_y-2 ; jj+=2) {
                for (kk=(N_z-D_absorb) ; kk<N_z-2 ; kk+=2) {      //NZ-d_absorb-2 ???
                    damp = ((double)kk-((double)N_z-(double)D_absorb))/(double)D_absorb;
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
                for (ii=2 ; ii<D_absorb-2 ; ii+=2) {
                    damp = ((double)ii-(double)D_absorb)/(double)D_absorb;
                    damp = ABSORBER_DAMPING(eco,damp);

                    EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                    EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                    EB_WAVE[ii  ][jj  ][kk+1] *= damp;
                }
            }
        }
    }
    if ( strstr(absorber,"x2") ) {
    // x2 absorber: x=d_absorb...NX
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (jj=2 ; jj<N_y-2 ; jj+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {  
                for (ii=(N_x-D_absorb) ; ii<N_x-2 ; ii+=2) {    //NX-d_absorb-2 ???
                    damp = ((double)ii-((double)N_x-(double)D_absorb))/(double)D_absorb;
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
                for (jj=2 ; jj<D_absorb-2 ; jj+=2) {
                    damp = ((double)jj-(double)D_absorb)/(double)D_absorb;
                    damp = ABSORBER_DAMPING(eco,damp);

                    EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                    EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                    EB_WAVE[ii  ][jj  ][kk+1] *= damp;
                }
            }
        }
    }
    // y2 absorber: y=d_absorb...NY
    if ( strstr(absorber,"y2") ) {
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {
                for (jj=(N_y-D_absorb) ; jj<N_y-2 ; jj+=2) {  //NY-d_absorb-2 ???
                    damp = ((double)jj-((double)N_y-(double)D_absorb))/(double)D_absorb;
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
                               double EB_WAVE[NX][NY][NZ] ) {
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
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (jj=2 ; jj<NY-2 ; jj+=2) {
            for (kk=2 ; kk<NZ-2 ; kk+=2) {
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
                           double EB_WAVE[NX][NY][NZ], 
                           double E_old[8][NY][NZ] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        jj, kk, 
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj<NY-2 ; jj+=2) {
        for (kk=2 ; kk<NZ-2 ; kk+=2) {
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

            // store values at x=NX-1 and x=NX-2
            // Ex: odd-even-even
            E_old[4+1][jj  ][kk  ]  = EB_WAVE[NX-4-offset+1][jj  ][kk  ];
            E_old[6+1][jj  ][kk  ]  = EB_WAVE[NX-2-offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_old[4  ][jj+1][kk  ]  = EB_WAVE[NX-4-offset  ][jj+1][kk  ];
            E_old[6  ][jj+1][kk  ]  = EB_WAVE[NX-2-offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_old[4  ][jj  ][kk+1]  = EB_WAVE[NX-4-offset  ][jj  ][kk+1];
            E_old[6  ][jj  ][kk+1]  = EB_WAVE[NX-2-offset  ][jj  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldE_ydir( gridConfiguration *gridCfg, 
                           double EB_WAVE[NX][NY][NZ], 
                           double E_old[NX][8][NZ] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, kk,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (kk=2 ; kk<NZ-2 ; kk+=2) {
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

            // store values at x=NX-1 and x=NX-2
            // Ex: odd-even-even
            E_old[ii+1][4  ][kk  ]  = EB_WAVE[ii+1][NY-4-offset  ][kk  ];
            E_old[ii+1][6  ][kk  ]  = EB_WAVE[ii+1][NY-2-offset  ][kk  ];
            // Ey: even-odd-even
            E_old[ii  ][4+1][kk  ]  = EB_WAVE[ii  ][NY-4-offset+1][kk  ];
            E_old[ii  ][6+1][kk  ]  = EB_WAVE[ii  ][NY-2-offset+1][kk  ];
            // Ez: even-even-odd
            E_old[ii  ][4  ][kk+1]  = EB_WAVE[ii  ][NY-4-offset  ][kk+1];
            E_old[ii  ][6  ][kk+1]  = EB_WAVE[ii  ][NY-2-offset  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldE_zdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[NX][NY][NZ], 
                           double E_old[NX][NY][8] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (jj=2 ; jj<NY-2 ; jj+=2) {
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

            // store values at z=NZ-1 and z=NZ-2
            // Ex: odd-even-even
            E_old[ii+1][jj  ][4  ]  = EB_WAVE[ii+1][jj  ][NZ-4-offset  ];
            E_old[ii+1][jj  ][6  ]  = EB_WAVE[ii+1][jj  ][NZ-2-offset  ];
            // Ey: even-odd-even
            E_old[ii  ][jj+1][4  ]  = EB_WAVE[ii  ][jj+1][NZ-4-offset  ];
            E_old[ii  ][jj+1][6  ]  = EB_WAVE[ii  ][jj+1][NZ-2-offset  ];
            // Ez: even-even-odd
            E_old[ii  ][jj  ][4+1]  = EB_WAVE[ii  ][jj  ][NZ-4-offset+1];
            E_old[ii  ][jj  ][6+1]  = EB_WAVE[ii  ][jj  ][NZ-2-offset+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldEref_xdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[NX][NY][NZ_REF], 
                              double E_old[8][NY][NZ_REF] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        jj, kk, 
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj<NY-2 ; jj+=2) {
        for (kk=2 ; kk<NZ_REF-2 ; kk+=2) {
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

            // store values at x=NX-1 and x=NX-2
            // Ex: odd-even-even
            E_old[4+1][jj  ][kk  ]  = EB_WAVE_ref[NX-4-offset+1][jj  ][kk  ];
            E_old[6+1][jj  ][kk  ]  = EB_WAVE_ref[NX-2-offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_old[4  ][jj+1][kk  ]  = EB_WAVE_ref[NX-4-offset  ][jj+1][kk  ];
            E_old[6  ][jj+1][kk  ]  = EB_WAVE_ref[NX-2-offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_old[4  ][jj  ][kk+1]  = EB_WAVE_ref[NX-4-offset  ][jj  ][kk+1];
            E_old[6  ][jj  ][kk+1]  = EB_WAVE_ref[NX-2-offset  ][jj  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldEref_ydir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[NX][NY][NZ_REF], 
                              double E_old[NX][8][NZ_REF] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, kk,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (kk=2 ; kk<NZ_REF-2 ; kk+=2) {
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

            // store values at x=NX-1 and x=NX-2
            // Ex: odd-even-even
            E_old[ii+1][4  ][kk  ]  = EB_WAVE_ref[ii+1][NY-4-offset  ][kk  ];
            E_old[ii+1][6  ][kk  ]  = EB_WAVE_ref[ii+1][NY-2-offset  ][kk  ];
            // Ey: even-odd-even
            E_old[ii  ][4+1][kk  ]  = EB_WAVE_ref[ii  ][NY-4-offset+1][kk  ];
            E_old[ii  ][6+1][kk  ]  = EB_WAVE_ref[ii  ][NY-2-offset+1][kk  ];
            // Ez: even-even-odd
            E_old[ii  ][4  ][kk+1]  = EB_WAVE_ref[ii  ][NY-4-offset  ][kk+1];
            E_old[ii  ][6  ][kk+1]  = EB_WAVE_ref[ii  ][NY-2-offset  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldEref_zdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[NX][NY][NZ_REF], 
                              double E_old[NX][NY][8] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (jj=2 ; jj<NY-2 ; jj+=2) {
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

            // store values at z=NZ-1 and z=NZ-2
            // Ex: odd-even-even
            E_old[ii+1][jj  ][4  ]  = EB_WAVE_ref[ii+1][jj  ][NZ_REF-4-offset  ];
            E_old[ii+1][jj  ][6  ]  = EB_WAVE_ref[ii+1][jj  ][NZ_REF-2-offset  ];
            // Ey: even-odd-even
            E_old[ii  ][jj+1][4  ]  = EB_WAVE_ref[ii  ][jj+1][NZ_REF-4-offset  ];
            E_old[ii  ][jj+1][6  ]  = EB_WAVE_ref[ii  ][jj+1][NZ_REF-2-offset  ];
            // Ez: even-even-odd
            E_old[ii  ][jj  ][4+1]  = EB_WAVE_ref[ii  ][jj  ][NZ_REF-4-offset+1];
            E_old[ii  ][jj  ][6+1]  = EB_WAVE_ref[ii  ][jj  ][NZ_REF-2-offset+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_1st( gridConfiguration *gridCfg, 
                 char absorber[],
                 double EB_WAVE[NX][NY][NZ], 
                 double E_old_xdir[8][NY][NZ], 
                 double E_old_ydir[NX][8][NZ], 
                 double E_old_zdir[NX][NY][8] ) {
//{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk,
        offset;             // refers to EB_WAVE only

    double
        cnst;

    cnst    = (DT-dx)/(DT+dx);
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
        for (jj=2 ; jj<NY-2 ; jj+=2) {
            for (kk=2 ; kk<NZ-2 ; kk+=2) {
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
        for (jj=2 ; jj<NY-2 ; jj+=2) {
            for (kk=2 ; kk<NZ-2 ; kk+=2) {
                // absorber at x=NX grid boundary
                // Ex: odd-even-even
                EB_WAVE[NX-2-offset+1][jj  ][kk  ]    = E_old_xdir[4+1][jj  ][kk  ]
                    + cnst * (    EB_WAVE[NX-4-offset+1][jj  ][kk  ] 
                              -E_old_xdir[6+1                   ][jj  ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[NX-2-offset  ][jj+1][kk  ]    = E_old_xdir[4  ][jj+1][kk  ]
                    + cnst * (    EB_WAVE[NX-4-offset  ][jj+1][kk  ] 
                              -E_old_xdir[6                     ][jj+1][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[NX-2-offset  ][jj  ][kk+1]    = E_old_xdir[4  ][jj  ][kk+1]
                    + cnst * (    EB_WAVE[NX-4-offset  ][jj  ][kk+1] 
                              -E_old_xdir[6                     ][jj  ][kk+1] );
            }
        }
    }

    // absorber into y-direction
    if ( strstr(absorber,"y1") ) {
        //printf("abs_Mur_1st_v2: y1\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
        for (ii=2 ; ii<NX-2 ; ii+=2) {
            for (kk=2 ; kk<NZ-2 ; kk+=2) {
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
        for (ii=2 ; ii<NX-2 ; ii+=2) {
            for (kk=2 ; kk<NZ-2 ; kk+=2) {
                // absorber at y=NY grid boundary
                // Ex: odd-even-even
                EB_WAVE[ii+1][NY-2-offset  ][kk  ] = E_old_ydir[ii+1][4  ][kk  ]
                    + cnst * (    EB_WAVE[ii+1][NY-4-offset  ][kk  ]
                              -E_old_ydir[ii+1][6                     ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[ii  ][NY-2-offset+1][kk  ] = E_old_ydir[ii  ][4+1][kk  ]
                    + cnst * (    EB_WAVE[ii  ][NY-4-offset+1][kk  ]
                              -E_old_ydir[ii  ][6+1                   ][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[ii  ][NY-2-offset  ][kk+1] = E_old_ydir[ii  ][4  ][kk+1]
                    + cnst * (    EB_WAVE[ii  ][NY-4-offset  ][kk+1]
                              -E_old_ydir[ii  ][6                     ][kk+1] );
            }
        }
    }

    // absorber into z-direction
    if ( strstr(absorber,"z1") ) {
        //printf("abs_Mur_1st_v2: z1\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
        for (ii=2 ; ii<NX-2 ; ii+=2) {
            for (jj=2 ; jj<NY-2 ; jj+=2) {
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
        for (ii=2 ; ii<NX-2 ; ii+=2) {
            for (jj=2 ; jj<NY-2 ; jj+=2) {
                // absorber at z=NZ grid boundary
                // Ex: odd-even-even
                EB_WAVE[ii+1][jj  ][NZ-2-offset  ]    = E_old_zdir[ii+1][jj  ][4  ]
                    + cnst * (    EB_WAVE[ii+1][jj  ][NZ-4-offset  ]
                              -E_old_zdir[ii+1][jj  ][6  ]     );
                // Ey: even-odd-even
                EB_WAVE[ii  ][jj+1][NZ-2-offset  ]    = E_old_zdir[ii  ][jj+1][4  ]
                    + cnst * (    EB_WAVE[ii  ][jj+1][NZ-4-offset  ]
                              -E_old_zdir[ii  ][jj+1][6  ]     );
                // Ez: even-even-odd
                EB_WAVE[ii  ][jj  ][NZ-2-offset+1]    = E_old_zdir[ii  ][jj  ][4+1]
                    + cnst * (    EB_WAVE[ii  ][jj  ][NZ-4-offset+1]
                              -E_old_zdir[ii  ][jj  ][6+1]     );
            }
        }
    }

    return EXIT_SUCCESS;

} //}}}


int abc_Mur_1st_ref( gridConfiguration *gridCfg,
                     double EB_WAVE[NX][NY][NZ_REF], 
                     double E_old_xdir[8][NY][NZ_REF], 
                     double E_old_ydir[NX][8][NZ_REF], 
                     double E_old_zdir[NX][NY][8] ) {
//{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk,
        offset;             // refers to EB_WAVE only

    double
        cnst;

    cnst    = (DT-dx)/(DT+dx);
    offset  = 2;

    // absorber into x-direction
#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj<NY-2 ; jj+=2) {
        for (kk=2 ; kk<NZ_REF-2 ; kk+=2) {
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
            // absorber at x=NX grid boundary
            // Ex: odd-even-even
            EB_WAVE[NX-2-offset+1][jj  ][kk  ]    = E_old_xdir[4+1][jj  ][kk  ]
                + cnst * (    EB_WAVE[NX-4-offset+1][jj  ][kk  ] 
                          -E_old_xdir[6+1                   ][jj  ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[NX-2-offset  ][jj+1][kk  ]    = E_old_xdir[4  ][jj+1][kk  ]
                + cnst * (    EB_WAVE[NX-4-offset  ][jj+1][kk  ] 
                          -E_old_xdir[6                     ][jj+1][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[NX-2-offset  ][jj  ][kk+1]    = E_old_xdir[4  ][jj  ][kk+1]
                + cnst * (    EB_WAVE[NX-4-offset  ][jj  ][kk+1] 
                          -E_old_xdir[6                     ][jj  ][kk+1] );
        }
    }

    // absorber into y-direction
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (kk=2 ; kk<NZ_REF-2 ; kk+=2) {
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
            // absorber at y=NY grid boundary
            // Ex: odd-even-even
            EB_WAVE[ii+1][NY-2-offset  ][kk  ] = E_old_ydir[ii+1][4  ][kk  ]
                + cnst * (    EB_WAVE[ii+1][NY-4-offset  ][kk  ]
                          -E_old_ydir[ii+1][6                     ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[ii  ][NY-2-offset+1][kk  ] = E_old_ydir[ii  ][4+1][kk  ]
                + cnst * (    EB_WAVE[ii  ][NY-4-offset+1][kk  ]
                          -E_old_ydir[ii  ][6+1                   ][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[ii  ][NY-2-offset  ][kk+1] = E_old_ydir[ii  ][4  ][kk+1]
                + cnst * (    EB_WAVE[ii  ][NY-4-offset  ][kk+1]
                          -E_old_ydir[ii  ][6                     ][kk+1] );
        }
    }

    // absorber into z-direction
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (jj=2 ; jj<NY-2 ; jj+=2) {
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
            // absorber at z=NZ grid boundary
            // Ex: odd-even-even
            EB_WAVE[ii+1][jj  ][NZ_REF-2-offset  ]    = E_old_zdir[ii+1][jj  ][4  ]
                + cnst * (    EB_WAVE[ii+1][jj  ][NZ_REF-4-offset  ]
                          -E_old_zdir[ii+1][jj  ][6  ]     );
            // Ey: even-odd-even
            EB_WAVE[ii  ][jj+1][NZ_REF-2-offset  ]    = E_old_zdir[ii  ][jj+1][4  ]
                + cnst * (    EB_WAVE[ii  ][jj+1][NZ_REF-4-offset  ]
                          -E_old_zdir[ii  ][jj+1][6  ]     );
            // Ez: even-even-odd
            EB_WAVE[ii  ][jj  ][NZ_REF-2-offset+1]    = E_old_zdir[ii  ][jj  ][4+1]
                + cnst * (    EB_WAVE[ii  ][jj  ][NZ_REF-4-offset+1]
                          -E_old_zdir[ii  ][jj  ][6+1]     );
        }
    }

    return EXIT_SUCCESS;

} //}}}


int set2zero_1D( size_t N_x, double arr_1D[N_x] ){
//{{{

    size_t
        ii;

#pragma omp parallel for default(shared) private(ii)
    for (ii=0 ; ii<N_x ; ++ii) {
        arr_1D[ii] = .0;
    }

    return EXIT_SUCCESS;
} //}}}


int set2zero_3D( size_t N_x, size_t N_y, size_t N_z, double arr_3D[N_x][N_y][N_z] ){
//{{{

    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii<N_x ; ++ii) {
        for (jj=0 ; jj<N_y ; ++jj) {
            for (kk=0 ; kk<N_z ; ++kk) {
                arr_3D[ii][jj][kk]  = .0;
            }
        }
    }

    return EXIT_SUCCESS;
} //}}}
