#include "focal.h"

void advance_fields(    gridConfiguration *gridCfg, 
                        double EB_WAVE[NX][NY][NZ], 
                        double EB_WAVE_ref[NX][NY][NZ_REF],
                        double J_B0[NX][NY][NZ],
                        double n_e[NX/2][NY/2][NZ/2]){
    //{{{
    
    if( BOUNDARY != 3 ){

        //Advance wave-plasma current
        advance_J(      gridCfg, EB_WAVE, J_B0, n_e );

        // advance wave magnetic field
        advance_B(      gridCfg, EB_WAVE );
        advance_B_ref(  gridCfg, EB_WAVE_ref );

        // advance wave electric field
        advance_E(      gridCfg, EB_WAVE, J_B0 );
        advance_E_ref(  gridCfg, EB_WAVE_ref );

    }else if( BOUNDARY == 3 ){

        //Advance wave-plasma current
        advance_J_UPML(      gridCfg, EB_WAVE, J_B0, n_e );

        //Advance wave magnetic field
        advance_B_UPML(      gridCfg, EB_WAVE);
        advance_Bref_UPML(   gridCfg, EB_WAVE_ref );

        //Advance wave electric field
        advance_E_UPML(      gridCfg, EB_WAVE, J_B0 );
        advance_Eref_UPML(   gridCfg, EB_WAVE_ref );

    }
}//}}}


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
                EB_WAVE[ii  ][jj+1][kk+1]   += -1.*DT/DX*(
                        +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                        );
                // -dBy/DT = dEx/dz - dEz/DX
                EB_WAVE[ii+1][jj  ][kk+1] += -1.*DT/DX*(
                        +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                        -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                        );
                // -dBz/DT = dEy/dx - dEx/dy
                EB_WAVE[ii+1][jj+1][kk  ] += -1.*DT/DX*(
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
                EB_WAVE[ii  ][jj+1][kk+1]   += -1.*DT/DX*(
                        +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                        );
                // -dBy/DT = dEx/dz - dEz/dx
                EB_WAVE[ii+1][jj  ][kk+1] += -1.*DT/DX*(
                        +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                        -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                        );
                // -dBz/DT = dEy/dx - dEx/dy
                EB_WAVE[ii+1][jj+1][kk  ] += -1.*DT/DX*(
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
                EB_WAVE[ii+1][jj  ][kk  ] += DT/DX*(
                        +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                        -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                        ) - DT*J_B0[ii+1][jj  ][kk  ];
                // dEy/DT = (dBx/dz - dBz/dx)
                EB_WAVE[ii  ][jj+1][kk  ] += DT/DX*(
                        +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                        -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                        ) - DT*J_B0[ii  ][jj+1][kk  ];
                // dEz/dt = (dBy/dx - dBx/dy)
                EB_WAVE[ii  ][jj  ][kk+1] += DT/DX*(
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
                EB_WAVE[ii+1][jj  ][kk  ] += DT/DX*(
                        +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                        -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                        );
                // dEy/dt = (dBx/dz - dBz/dx)
                EB_WAVE[ii  ][jj+1][kk  ] += DT/DX*(
                        +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                        -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                        );
                // dEz/dt = (dBy/dx - dBx/dy)
                EB_WAVE[ii  ][jj  ][kk+1] += DT/DX*(
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
    ne_dist     = round( PERIOD/1 );

    printf( "Starting to set density in absorber to %f with a smooth transition...\n", ne_absorb);

    x0          = (double)D_ABSORB + ne_dist;
    x1          = (double)NX - (D_ABSORB + ne_dist);
    y0          = (double)D_ABSORB + ne_dist;
    y1          = (double)NY - (D_ABSORB + ne_dist);
    z0          = (double)D_ABSORB + ne_dist;
    z1          = (double)NZ - (D_ABSORB + ne_dist);

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

    printf( "...done setting density in absorber to %f\n", ne_absorb);

    return EXIT_SUCCESS;
} //}}}


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

/*UPML functions*/
int advance_J_UPML(  gridConfiguration *gridCfg, 
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
    for (ii= D_ABSORB ; ii<=NX-2-D_ABSORB ; ii+=2) {
        for (jj= D_ABSORB ; jj<=NY-2-D_ABSORB ; jj+=2) {
            for (kk= D_ABSORB ; kk<=NZ-2-D_ABSORB ; kk+=2) {
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


int advance_B_UPML( gridConfiguration *gridCfg, 
                    double EB_WAVE[NX][NY][NZ] ) {
//{{{
    // B_new = B_old - nabla x E

    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=D_ABSORB ; ii<=NX-2-D_ABSORB ; ii+=2) {
        for (jj=D_ABSORB ; jj<=NY-2-D_ABSORB ; jj+=2) {
            for (kk=D_ABSORB ; kk<=NZ-2-D_ABSORB ; kk+=2) {

                // -dBx/DT = dEz/dy - dEy/dz
                EB_WAVE[ii  ][jj+1][kk+1]   += -1.*DT/DX*(
                        +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                        );
                // -dBy/DT = dEx/dz - dEz/DX
                EB_WAVE[ii+1][jj  ][kk+1] += -1.*DT/DX*(
                        +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                        -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                        );
                // -dBz/DT = dEy/dx - dEx/dy
                EB_WAVE[ii+1][jj+1][kk  ] += -1.*DT/DX*(
                        +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                        -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                        );

            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_Bref_UPML(  gridConfiguration *gridCfg, 
                        double EB_WAVE[NX][NY][NZ_REF] ) {
//{{{
    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=D_ABSORB ; ii<=NX-2-D_ABSORB ; ii+=2) {
        for (jj=D_ABSORB ; jj<=NY-2-D_ABSORB ; jj+=2) {
            for (kk=D_ABSORB ; kk<=NZ_REF-2-D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                EB_WAVE[ii  ][jj+1][kk+1]   += -1.*DT/DX*(
                        +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                        );
                // -dBy/DT = dEx/dz - dEz/dx
                EB_WAVE[ii+1][jj  ][kk+1] += -1.*DT/DX*(
                        +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                        -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                        );
                // -dBz/DT = dEy/dx - dEx/dy
                EB_WAVE[ii+1][jj+1][kk  ] += -1.*DT/DX*(
                        +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                        -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                        );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_E_UPML(  gridConfiguration *gridCfg, 
                    double EB_WAVE[NX][NY][NZ], 
                    double J_B0[NX][NY][NZ] ) {
//{{{
    // E_new = E_old + c^2*nablaxB - 1/epsilon_0*J

    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=D_ABSORB ; ii<=NX-2-D_ABSORB ; ii+=2) {
        for (jj=D_ABSORB ; jj<=NY-2-D_ABSORB ; jj+=2) {
            for (kk=D_ABSORB ; kk<=NZ-2-D_ABSORB ; kk+=2) {
                // dEx/DT = (dBz/dy - dBy/dz)
                EB_WAVE[ii+1][jj  ][kk  ] += DT/DX*(
                        +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                        -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                        ) - DT*J_B0[ii+1][jj  ][kk  ];
                // dEy/DT = (dBx/dz - dBz/dx)
                EB_WAVE[ii  ][jj+1][kk  ] += DT/DX*(
                        +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                        -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                        ) - DT*J_B0[ii  ][jj+1][kk  ];
                // dEz/dt = (dBy/dx - dBx/dy)
                EB_WAVE[ii  ][jj  ][kk+1] += DT/DX*(
                        +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                        ) - DT*J_B0[ii  ][jj  ][kk+1];

            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_Eref_UPML(   gridConfiguration *gridCfg, 
                        double EB_WAVE[NX][NY][NZ_REF] ) { 
//{{{
    // same as advance_E but for reference fields (directional coupler)
    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=D_ABSORB ; ii<=NX-2-D_ABSORB ; ii+=2) {
        for (jj=D_ABSORB ; jj<=NY-2-D_ABSORB ; jj+=2) {
            for (kk=D_ABSORB ; kk<=NZ_REF-2-D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                EB_WAVE[ii+1][jj  ][kk  ] += DT/DX*(
                        +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                        -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                        );
                // dEy/dt = (dBx/dz - dBz/dx)
                EB_WAVE[ii  ][jj+1][kk  ] += DT/DX*(
                        +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                        -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                        );
                // dEz/dt = (dBy/dx - dBx/dy)
                EB_WAVE[ii  ][jj  ][kk+1] += DT/DX*(
                        +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                        );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}
