#include "power_calc.h"

static double **timetraces = NULL;

//Initialization of the power related values
int init_powerValues(   gridConfiguration *gridCfg,
                        powerValues *powerVal ){

    // the arrays are initialized with calloc() and thus don't require zeroing
    timetraces = allocate2DArray( (T_END/(int)PERIOD), 8 );

    printf( "starting to set all variables to 0...\n" );
    power_abs_x1    = .0;
    power_abs_x2    = .0;
    power_abs_y1    = .0;
    power_abs_y2    = .0;
    power_abs_z1    = .0;
    power_abs_z2    = .0;
    power_abs_ref   = 1e-7;
    poynt_x1       = .0;
    poynt_x2       = .0;
    poynt_y1       = .0;
    poynt_y2       = .0;
    poynt_z1       = .0;
    poynt_z1_ref   = .0;
    poynt_z2       = .0;
    printf( "...done setting all variables to 0\n" );

    pwr_dect    = D_ABSORB;

    return EXIT_SUCCESS;
}

//function called in MAIN to compute power
void compute_power( gridConfiguration *gridCfg, 
                    beamAntennaConfiguration *beamCfg,
                    powerValues *powerVal,
                    int t_int, 
                    double EB_WAVE[NX][NY][NZ], 
                    double EB_WAVE_ref[NX][NY][NZ_REF] ){

    calculate_power( gridCfg, powerVal, t_int, EB_WAVE, EB_WAVE_ref );
    power_toTimetraces( gridCfg, beamCfg, powerVal, t_int );

}

//calls power calculation
int calculate_power(    gridConfiguration *gridCfg, 
                        powerValues *powerVal,
                        int t_int, 
                        double EB_WAVE[NX][NY][NZ], 
                        double EB_WAVE_ref[NX][NY][NZ_REF] ){

    // IQ detector for power detection
    if ( t_int >= 20*PERIOD ) {
        // z1-plane and z2-plane
        poynt_z1_ref    = calc_poynt_4( gridCfg, powerVal, "ref_z1", EB_WAVE, EB_WAVE_ref );
        poynt_z1        = calc_poynt_4( gridCfg, powerVal, "z1",     EB_WAVE, EB_WAVE_ref );
        poynt_z2        = calc_poynt_4( gridCfg, powerVal, "z2",     EB_WAVE, EB_WAVE_ref );
        // x1-plane and x2-plane
        poynt_x1        = calc_poynt_4( gridCfg, powerVal, "x1", EB_WAVE, EB_WAVE_ref );
        poynt_x2        = calc_poynt_4( gridCfg, powerVal, "x2", EB_WAVE, EB_WAVE_ref );
        // y1-plane and y2-plane
        poynt_y1        = calc_poynt_4( gridCfg, powerVal, "y1", EB_WAVE, EB_WAVE_ref );
        poynt_y2        = calc_poynt_4( gridCfg, powerVal, "y2", EB_WAVE, EB_WAVE_ref );
            
        //printf( "t = %d, power_abs_ref = %13.5e, power_abs_z1 = %13.5e, power_abs_z2 = %13.5e, poynt_z1 = %13.5e, poynt_z2 = %13.5e\n",
        //        t_int, power_abs_ref, power_abs_z1, power_abs_z2, poynt_z1, poynt_z2 );

        power_abs_ref   = .99*power_abs_ref + .01*poynt_z1_ref;
        power_abs_z1    = .99*power_abs_z1  + .01*poynt_z1;
        power_abs_z2    = .99*power_abs_z2  + .01*poynt_z2;
        power_abs_x1    = .99*power_abs_x1  + .01*poynt_x1;
        power_abs_x2    = .99*power_abs_x2  + .01*poynt_x2;
        power_abs_y1    = .99*power_abs_y1  + .01*poynt_y1;
        power_abs_y2    = .99*power_abs_y2  + .01*poynt_y2;
    }

    return EXIT_SUCCESS;
}

//Store power values in Timetraces
int power_toTimetraces( gridConfiguration *gridCfg, 
                        beamAntennaConfiguration *beamCfg, 
                        powerValues *powerVal,
                        int t_int ){

    if ( (t_int % (int)(PERIOD)) == 4 )  {
        //print to console
        printf( "status: number of oscillation periods: %d (t_int= %d) \n",T_WAVE,t_int);
        printf( "        Poynting-power: z1 = %13.6e, z2 = %13.6e, x1 = %13.6e, x2 = %13.6e, y1 = %13.6e, y2 = %13.6e, (z1+z2+x1+x2+y1+y2)/z1_ref = %13.6e %%\n",
                power_abs_z1/power_abs_ref, 
                power_abs_z2/power_abs_ref,
                power_abs_x1/power_abs_ref, 
                power_abs_x2/power_abs_ref,
                power_abs_y1/power_abs_ref, 
                power_abs_y2/power_abs_ref,
                (power_abs_x1+power_abs_x2 + power_abs_y1+power_abs_y2 + power_abs_z1+power_abs_z2)/power_abs_ref * 100.
            );
        
        //Save to time traces
        timetraces[T_WAVE][0]   = (double)t_int;
        timetraces[T_WAVE][1]   = (double)T_WAVE;
        timetraces[T_WAVE][2]   = power_abs_z1/power_abs_ref;
        timetraces[T_WAVE][3]   = power_abs_z2/power_abs_ref;
        timetraces[T_WAVE][4]   = power_abs_x1/power_abs_ref;
        timetraces[T_WAVE][5]   = power_abs_x2/power_abs_ref;
        timetraces[T_WAVE][6]   = power_abs_y1/power_abs_ref;
        timetraces[T_WAVE][7]   = power_abs_y2/power_abs_ref;

    }

    return EXIT_SUCCESS;
}

//functions to compute Pointyng vector
double calc_poynt_4( gridConfiguration *gridCfg, 
                     powerValues *powerVal, 
                     char absorber[],
                     double EB_WAVE[NX][NY][NZ], 
                     double EB_WAVE_ref[NX][NY][NZ_REF] ) {
//{{{

    size_t
        ii, jj, kk;
    double
        poynt;

    poynt   = .0;

    // P = E x H
    // Px = Ey*Hz - Ez*Hy
    // Py = Ez*Hx - Ex*Hz
    // Pz = Ex*Hy - Ey*Hx
    // Bx: even-odd-odd
    // By: odd-even-odd
    // Bz: odd-odd-even
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd
    
    if ( strcmp(absorber,"ref_z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(NX-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<=(NY-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE_ref[ii+1][jj  ][pwr_dect  ]
                          *EB_WAVE_ref[ii+1][jj  ][pwr_dect+1]
                          -EB_WAVE_ref[ii  ][jj+1][pwr_dect  ]
                          *EB_WAVE_ref[ii  ][jj+1][pwr_dect+1] );
            }
        }
    } else if ( strcmp(absorber,"z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(NX-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<=(NY-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( ( EB_WAVE[ii+1][jj  ][pwr_dect  ] - EB_WAVE_ref[ii+1][jj  ][pwr_dect  ] )
                          *( EB_WAVE[ii+1][jj  ][pwr_dect+1] - EB_WAVE_ref[ii+1][jj  ][pwr_dect+1] )
                          -( EB_WAVE[ii  ][jj+1][pwr_dect  ] - EB_WAVE_ref[ii  ][jj+1][pwr_dect  ] ) 
                          *( EB_WAVE[ii  ][jj+1][pwr_dect+1] - EB_WAVE_ref[ii  ][jj+1][pwr_dect+1] ) );
            }
        }
    } else if ( strcmp(absorber,"z2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(NX-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<=(NY-pwr_dect-2) ; jj+=2) {
                // z2-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE[ii+1][jj  ][NZ-pwr_dect  ]
                          *EB_WAVE[ii+1][jj  ][NZ-pwr_dect-1]
                          -EB_WAVE[ii  ][jj+1][NZ-pwr_dect  ]
                          *EB_WAVE[ii  ][jj+1][NZ-pwr_dect-1] );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(NY-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(NZ-pwr_dect-2) ; kk+=2) {
                // x1-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[pwr_dect  ][jj+1][kk  ]
                          *EB_WAVE[pwr_dect+1][jj+1][kk  ]
                          -EB_WAVE[pwr_dect  ][jj  ][kk+1]
                          *EB_WAVE[pwr_dect+1][jj  ][kk+1] );
            }
        }
    } else if ( strcmp(absorber,"x2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(NY-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(NZ-pwr_dect-2) ; kk+=2) {
                // x2-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[NX-pwr_dect  ][jj+1][kk  ]
                          *EB_WAVE[NX-pwr_dect-1][jj+1][kk  ]
                          -EB_WAVE[NX-pwr_dect  ][jj  ][kk+1]
                          *EB_WAVE[NX-pwr_dect-1][jj  ][kk+1] );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(NX-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(NZ-pwr_dect-2) ; kk+=2) {
                // y1-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][pwr_dect  ][kk+1]
                          *EB_WAVE[ii  ][pwr_dect+1][kk+1]
                          -EB_WAVE[ii+1][pwr_dect  ][kk  ]
                          *EB_WAVE[ii+1][pwr_dect+1][kk  ] );
            }
        }
    } else if ( strcmp(absorber,"y2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(NX-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(NZ-pwr_dect-2) ; kk+=2) {
                // y2-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][NY-pwr_dect  ][kk+1]
                          *EB_WAVE[ii  ][NY-pwr_dect-1][kk+1]
                          -EB_WAVE[ii+1][NY-pwr_dect  ][kk  ]
                          *EB_WAVE[ii+1][NY-pwr_dect-1][kk  ] );
            }
        }
    }
    
    return fabs(poynt);
} //}}}


double calc_poynt_5( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     powerValues *powerVal, 
                     char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] ) {
//{{{
// same as calc_poynt_4, but with spatial averaging of B-component

    size_t
        ii, jj, kk;
    double
        poynt;

    poynt   = .0;

    // P = E x H
    // Px = Ey*Hz - Ez*Hy
    // Py = Ez*Hx - Ex*Hz
    // Pz = Ex*Hy - Ey*Hx
    // Bx: even-odd-odd
    // By: odd-even-odd
    // Bz: odd-odd-even
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd
    
    if ( strcmp(absorber,"ref_z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE_ref[ii+1][jj  ][pwr_dect  ]
                          *( EB_WAVE_ref[ii-1][jj  ][pwr_dect+1]
                            +EB_WAVE_ref[ii+3][jj  ][pwr_dect+1]
                            +EB_WAVE_ref[ii+1][jj-2][pwr_dect+1]
                            +EB_WAVE_ref[ii+1][jj+2][pwr_dect+1] )*.25
                          -EB_WAVE_ref[ii  ][jj+1][pwr_dect  ]
                          *( EB_WAVE_ref[ii-2][jj+1][pwr_dect+1]
                            +EB_WAVE_ref[ii+2][jj+1][pwr_dect+1]
                            +EB_WAVE_ref[ii  ][jj-1][pwr_dect+1]
                            +EB_WAVE_ref[ii  ][jj+3][pwr_dect+1] )*.25 );
            }
        }
    } else if ( strcmp(absorber,"z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( ( EB_WAVE[ii+1][jj  ][pwr_dect  ] - EB_WAVE_ref[ii+1][jj  ][pwr_dect  ] )
                          *( ( EB_WAVE[ii-1][jj  ][pwr_dect+1] 
                              +EB_WAVE[ii+3][jj  ][pwr_dect+1] 
                              +EB_WAVE[ii+1][jj-2][pwr_dect+1] 
                              +EB_WAVE[ii+1][jj+2][pwr_dect+1] )*.25
                            -( EB_WAVE_ref[ii-1][jj  ][pwr_dect+1]
                              +EB_WAVE_ref[ii+3][jj  ][pwr_dect+1]
                              +EB_WAVE_ref[ii+1][jj-2][pwr_dect+1]
                              +EB_WAVE_ref[ii+1][jj+2][pwr_dect+1] )*.25 )
                          -( EB_WAVE[ii  ][jj+1][pwr_dect  ] - EB_WAVE_ref[ii  ][jj+1][pwr_dect  ] ) 
                          *( ( EB_WAVE[ii-2][jj+1][pwr_dect+1] 
                              +EB_WAVE[ii+2][jj+1][pwr_dect+1] 
                              +EB_WAVE[ii  ][jj-1][pwr_dect+1] 
                              +EB_WAVE[ii  ][jj+3][pwr_dect+1] )*.25
                            -( EB_WAVE_ref[ii-2][jj+1][pwr_dect+1]
                              +EB_WAVE_ref[ii+2][jj+1][pwr_dect+1]
                              +EB_WAVE_ref[ii  ][jj-1][pwr_dect+1]
                              +EB_WAVE_ref[ii  ][jj+3][pwr_dect+1] )*.25 ) );
            }
        }
    } else if ( strcmp(absorber,"z2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z2-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE[ii+1][jj  ][N_z-pwr_dect  ]
                          *( EB_WAVE[ii-1][jj  ][N_z-pwr_dect-1]
                            +EB_WAVE[ii+3][jj  ][N_z-pwr_dect-1]
                            +EB_WAVE[ii+1][jj-2][N_z-pwr_dect-1]
                            +EB_WAVE[ii+1][jj+2][N_z-pwr_dect-1] )*.25
                          -EB_WAVE[ii  ][jj+1][N_z-pwr_dect  ]
                          *( EB_WAVE[ii-2][jj+1][N_z-pwr_dect-1]
                            +EB_WAVE[ii+2][jj+1][N_z-pwr_dect-1]
                            +EB_WAVE[ii  ][jj-1][N_z-pwr_dect-1]
                            +EB_WAVE[ii  ][jj+3][N_z-pwr_dect-1] )*.25 );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x1-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[pwr_dect  ][jj+1][kk  ]
                          *( EB_WAVE[pwr_dect+1][jj-1][kk  ]
                            +EB_WAVE[pwr_dect+1][jj+3][kk  ]
                            +EB_WAVE[pwr_dect+1][jj+1][kk-2]
                            +EB_WAVE[pwr_dect+1][jj+1][kk+2] )*.25
                          -EB_WAVE[pwr_dect  ][jj  ][kk+1]
                          *( EB_WAVE[pwr_dect+1][jj-2][kk+1]
                            +EB_WAVE[pwr_dect+1][jj+2][kk+1]
                            +EB_WAVE[pwr_dect+1][jj  ][kk-1]
                            +EB_WAVE[pwr_dect+1][jj  ][kk+3] )*.25 );
            }
        }
    } else if ( strcmp(absorber,"x2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x2-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[N_x-pwr_dect  ][jj+1][kk  ]
                          *( EB_WAVE[N_x-pwr_dect-1][jj-1][kk  ]
                            +EB_WAVE[N_x-pwr_dect-1][jj+3][kk  ]
                            +EB_WAVE[N_x-pwr_dect-1][jj+1][kk-2]
                            +EB_WAVE[N_x-pwr_dect-1][jj+1][kk+2] )*.25
                          -EB_WAVE[N_x-pwr_dect  ][jj  ][kk+1]
                          *( EB_WAVE[N_x-pwr_dect-1][jj-2][kk+1]
                            +EB_WAVE[N_x-pwr_dect-1][jj+2][kk+1]
                            +EB_WAVE[N_x-pwr_dect-1][jj  ][kk-1]
                            +EB_WAVE[N_x-pwr_dect-1][jj  ][kk+3] )*.25 );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y1-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][pwr_dect  ][kk+1]
                          *( EB_WAVE[ii-2][pwr_dect+1][kk+1]
                            +EB_WAVE[ii+2][pwr_dect+1][kk+1]
                            +EB_WAVE[ii  ][pwr_dect+1][kk-1]
                            +EB_WAVE[ii  ][pwr_dect+1][kk+3] )*.25
                          -EB_WAVE[ii+1][pwr_dect  ][kk  ]
                          *( EB_WAVE[ii-1][pwr_dect+1][kk  ]
                            +EB_WAVE[ii+3][pwr_dect+1][kk  ]
                            +EB_WAVE[ii+1][pwr_dect+1][kk-2]
                            +EB_WAVE[ii+1][pwr_dect+1][kk+2] )*.25 );
            }
        }
    } else if ( strcmp(absorber,"y2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y2-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][N_y-pwr_dect  ][kk+1]
                          *( EB_WAVE[ii-2][N_y-pwr_dect-1][kk+1]
                            +EB_WAVE[ii+2][N_y-pwr_dect-1][kk+1]
                            +EB_WAVE[ii  ][N_y-pwr_dect-1][kk-1]
                            +EB_WAVE[ii  ][N_y-pwr_dect-1][kk+3] )*.25
                          -EB_WAVE[ii+1][N_y-pwr_dect  ][kk  ]
                          *( EB_WAVE[ii-1][N_y-pwr_dect-1][kk  ]
                            +EB_WAVE[ii+3][N_y-pwr_dect-1][kk  ]
                            +EB_WAVE[ii+1][N_y-pwr_dect-1][kk-2]
                            +EB_WAVE[ii+1][N_y-pwr_dect-1][kk+2] )*.25 );
            }
        }
    }
    
    return fabs(poynt);
} //}}}


double calc_poynt_6( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     powerValues *powerVal, 
                     char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] ) {
//{{{
// same as calc_poynt_5, but full spatial averaging of B-component

    size_t
        ii, jj, kk;
    double
        poynt;

    poynt   = .0;

    // P = E x H
    // Px = Ey*Hz - Ez*Hy
    // Py = Ez*Hx - Ex*Hz
    // Pz = Ex*Hy - Ey*Hx
    // Bx: even-odd-odd
    // By: odd-even-odd
    // Bz: odd-odd-even
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd
    
    if ( strcmp(absorber,"ref_z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE_ref[ii+1][jj  ][pwr_dect  ]
                          *( EB_WAVE_ref[ii-1][jj  ][pwr_dect+1]
                            +EB_WAVE_ref[ii+3][jj  ][pwr_dect+1]
                            +EB_WAVE_ref[ii+1][jj-2][pwr_dect+1]
                            +EB_WAVE_ref[ii+1][jj+2][pwr_dect+1] 
                            +EB_WAVE_ref[ii+1][jj  ][pwr_dect-1]
                            +EB_WAVE_ref[ii+1][jj  ][pwr_dect+3] )/6.
                          -EB_WAVE_ref[ii  ][jj+1][pwr_dect  ]
                          *( EB_WAVE_ref[ii-2][jj+1][pwr_dect+1]
                            +EB_WAVE_ref[ii+2][jj+1][pwr_dect+1]
                            +EB_WAVE_ref[ii  ][jj-1][pwr_dect+1]
                            +EB_WAVE_ref[ii  ][jj+3][pwr_dect+1]
                            +EB_WAVE_ref[ii  ][jj+1][pwr_dect-1]
                            +EB_WAVE_ref[ii  ][jj+1][pwr_dect+3] )/6. );
            }
        }
    } else if ( strcmp(absorber,"z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( ( EB_WAVE[ii+1][jj  ][pwr_dect  ] - EB_WAVE_ref[ii+1][jj  ][pwr_dect  ] )
                          *( ( EB_WAVE[ii-1][jj  ][pwr_dect+1] 
                              +EB_WAVE[ii+3][jj  ][pwr_dect+1] 
                              +EB_WAVE[ii+1][jj-2][pwr_dect+1] 
                              +EB_WAVE[ii+1][jj+2][pwr_dect+1]
                              +EB_WAVE[ii+1][jj  ][pwr_dect-1] 
                              +EB_WAVE[ii+1][jj  ][pwr_dect+3] )/6.
                            -( EB_WAVE_ref[ii-1][jj  ][pwr_dect+1]
                              +EB_WAVE_ref[ii+3][jj  ][pwr_dect+1]
                              +EB_WAVE_ref[ii+1][jj-2][pwr_dect+1]
                              +EB_WAVE_ref[ii+1][jj+2][pwr_dect+1]
                              +EB_WAVE_ref[ii+1][jj  ][pwr_dect-1]
                              +EB_WAVE_ref[ii+1][jj  ][pwr_dect+3] )/6. )
                          -( EB_WAVE[ii  ][jj+1][pwr_dect  ] - EB_WAVE_ref[ii  ][jj+1][pwr_dect  ] ) 
                          *( ( EB_WAVE[ii-2][jj+1][pwr_dect+1] 
                              +EB_WAVE[ii+2][jj+1][pwr_dect+1] 
                              +EB_WAVE[ii  ][jj-1][pwr_dect+1] 
                              +EB_WAVE[ii  ][jj+3][pwr_dect+1]
                              +EB_WAVE[ii  ][jj+1][pwr_dect-1] 
                              +EB_WAVE[ii  ][jj+1][pwr_dect+3] )/6.
                            -( EB_WAVE_ref[ii-2][jj+1][pwr_dect+1]
                              +EB_WAVE_ref[ii+2][jj+1][pwr_dect+1]
                              +EB_WAVE_ref[ii  ][jj-1][pwr_dect+1]
                              +EB_WAVE_ref[ii  ][jj+3][pwr_dect+1]
                              +EB_WAVE_ref[ii  ][jj+1][pwr_dect-1]
                              +EB_WAVE_ref[ii  ][jj+1][pwr_dect+3] )/6. ) );
            }
        }
    } else if ( strcmp(absorber,"z2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z2-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE[ii+1][jj  ][N_z-pwr_dect  ]
                          *( EB_WAVE[ii-1][jj  ][N_z-pwr_dect-1]
                            +EB_WAVE[ii+3][jj  ][N_z-pwr_dect-1]
                            +EB_WAVE[ii+1][jj-2][N_z-pwr_dect-1]
                            +EB_WAVE[ii+1][jj+2][N_z-pwr_dect-1]
                            +EB_WAVE[ii+1][jj  ][N_z-pwr_dect-3]
                            +EB_WAVE[ii+1][jj  ][N_z-pwr_dect+1] )/6.
                          -EB_WAVE[ii  ][jj+1][N_z-pwr_dect  ]
                          *( EB_WAVE[ii-2][jj+1][N_z-pwr_dect-1]
                            +EB_WAVE[ii+2][jj+1][N_z-pwr_dect-1]
                            +EB_WAVE[ii  ][jj-1][N_z-pwr_dect-1]
                            +EB_WAVE[ii  ][jj+3][N_z-pwr_dect-1]
                            +EB_WAVE[ii  ][jj+1][N_z-pwr_dect-3]
                            +EB_WAVE[ii  ][jj+1][N_z-pwr_dect+1] )/6. );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x1-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[pwr_dect  ][jj+1][kk  ]
                          *( EB_WAVE[pwr_dect+1][jj-1][kk  ]
                            +EB_WAVE[pwr_dect+1][jj+3][kk  ]
                            +EB_WAVE[pwr_dect+1][jj+1][kk-2]
                            +EB_WAVE[pwr_dect+1][jj+1][kk+2]
                            +EB_WAVE[pwr_dect-1][jj+1][kk  ]
                            +EB_WAVE[pwr_dect+3][jj+1][kk  ] )/6.
                          -EB_WAVE[pwr_dect  ][jj  ][kk+1]
                          *( EB_WAVE[pwr_dect+1][jj-2][kk+1]
                            +EB_WAVE[pwr_dect+1][jj+2][kk+1]
                            +EB_WAVE[pwr_dect+1][jj  ][kk-1]
                            +EB_WAVE[pwr_dect+1][jj  ][kk+3]
                            +EB_WAVE[pwr_dect-1][jj  ][kk+1]
                            +EB_WAVE[pwr_dect+3][jj  ][kk+1] )/6. );
            }
        }
    } else if ( strcmp(absorber,"x2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x2-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[N_x-pwr_dect  ][jj+1][kk  ]
                          *( EB_WAVE[N_x-pwr_dect-1][jj-1][kk  ]
                            +EB_WAVE[N_x-pwr_dect-1][jj+3][kk  ]
                            +EB_WAVE[N_x-pwr_dect-1][jj+1][kk-2]
                            +EB_WAVE[N_x-pwr_dect-1][jj+1][kk+2]
                            +EB_WAVE[N_x-pwr_dect-3][jj+1][kk  ]
                            +EB_WAVE[N_x-pwr_dect+1][jj+1][kk  ] )/6.
                          -EB_WAVE[N_x-pwr_dect  ][jj  ][kk+1]
                          *( EB_WAVE[N_x-pwr_dect-1][jj-2][kk+1]
                            +EB_WAVE[N_x-pwr_dect-1][jj+2][kk+1]
                            +EB_WAVE[N_x-pwr_dect-1][jj  ][kk-1]
                            +EB_WAVE[N_x-pwr_dect-1][jj  ][kk+3]
                            +EB_WAVE[N_x-pwr_dect-3][jj  ][kk+1]
                            +EB_WAVE[N_x-pwr_dect+1][jj  ][kk+1] )/6. );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y1-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][pwr_dect  ][kk+1]
                          *( EB_WAVE[ii-2][pwr_dect+1][kk+1]
                            +EB_WAVE[ii+2][pwr_dect+1][kk+1]
                            +EB_WAVE[ii  ][pwr_dect+1][kk-1]
                            +EB_WAVE[ii  ][pwr_dect+1][kk+3]
                            +EB_WAVE[ii  ][pwr_dect-1][kk+1]
                            +EB_WAVE[ii  ][pwr_dect+3][kk+1] )/6.
                          -EB_WAVE[ii+1][pwr_dect  ][kk  ]
                          *( EB_WAVE[ii-1][pwr_dect+1][kk  ]
                            +EB_WAVE[ii+3][pwr_dect+1][kk  ]
                            +EB_WAVE[ii+1][pwr_dect+1][kk-2]
                            +EB_WAVE[ii+1][pwr_dect+1][kk+2]
                            +EB_WAVE[ii+1][pwr_dect-1][kk  ]
                            +EB_WAVE[ii+1][pwr_dect+3][kk  ] )/6. );
            }
        }
    } else if ( strcmp(absorber,"y2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y2-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][N_y-pwr_dect  ][kk+1]
                          *( EB_WAVE[ii-2][N_y-pwr_dect-1][kk+1]
                            +EB_WAVE[ii+2][N_y-pwr_dect-1][kk+1]
                            +EB_WAVE[ii  ][N_y-pwr_dect-1][kk-1]
                            +EB_WAVE[ii  ][N_y-pwr_dect-1][kk+3]
                            +EB_WAVE[ii  ][N_y-pwr_dect-3][kk+1]
                            +EB_WAVE[ii  ][N_y-pwr_dect+1][kk+1] )/6.
                          -EB_WAVE[ii+1][N_y-pwr_dect  ][kk  ]
                          *( EB_WAVE[ii-1][N_y-pwr_dect-1][kk  ]
                            +EB_WAVE[ii+3][N_y-pwr_dect-1][kk  ]
                            +EB_WAVE[ii+1][N_y-pwr_dect-1][kk-2]
                            +EB_WAVE[ii+1][N_y-pwr_dect-1][kk+2]
                            +EB_WAVE[ii+1][N_y-pwr_dect-3][kk  ]
                            +EB_WAVE[ii+1][N_y-pwr_dect+1][kk  ] )/6. );
            }
        }
    }
    
    return fabs(poynt);
} //}}}


double calc_poynt_7( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     powerValues *powerVal, 
                     char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] ) {
//{{{
// same as calc_poynt_4, but with symmetrizing of B-component

    size_t
        ii, jj, kk;
    double
        poynt;

    poynt   = .0;

    // P = E x H
    // Px = Ey*Hz - Ez*Hy
    // Py = Ez*Hx - Ex*Hz
    // Pz = Ex*Hy - Ey*Hx
    // Bx: even-odd-odd
    // By: odd-even-odd
    // Bz: odd-odd-even
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd
    
    if ( strcmp(absorber,"ref_z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE_ref[ii+1][jj  ][pwr_dect  ]
                          *( EB_WAVE_ref[ii+1][jj  ][pwr_dect+1]
                            +EB_WAVE_ref[ii+1][jj  ][pwr_dect-1] )*.5
                          -EB_WAVE_ref[ii  ][jj+1][pwr_dect  ]
                          *( EB_WAVE_ref[ii  ][jj+1][pwr_dect+1]
                            +EB_WAVE_ref[ii  ][jj+1][pwr_dect-1] )*.5 );
            }
        }
    } else if ( strcmp(absorber,"z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( ( EB_WAVE[ii+1][jj  ][pwr_dect  ] - EB_WAVE_ref[ii+1][jj  ][pwr_dect  ] )
                          *( ( EB_WAVE[ii+1][jj  ][pwr_dect+1] 
                              +EB_WAVE[ii+1][jj  ][pwr_dect-1] )*.5
                            -( EB_WAVE_ref[ii+1][jj  ][pwr_dect+1]
                              +EB_WAVE_ref[ii+1][jj  ][pwr_dect-1] )*.5 )
                          -( EB_WAVE[ii  ][jj+1][pwr_dect  ] - EB_WAVE_ref[ii  ][jj+1][pwr_dect  ] ) 
                          *( ( EB_WAVE[ii  ][jj+1][pwr_dect+1] 
                              +EB_WAVE[ii  ][jj+1][pwr_dect-1] )*.5
                            -( EB_WAVE_ref[ii  ][jj+1][pwr_dect+1]
                              +EB_WAVE_ref[ii  ][jj+1][pwr_dect-1] )*.5 ) );
            }
        }
    } else if ( strcmp(absorber,"z2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z2-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE[ii+1][jj  ][N_z-pwr_dect  ]
                          *( EB_WAVE[ii+1][jj  ][N_z-pwr_dect-1]
                            +EB_WAVE[ii+1][jj  ][N_z-pwr_dect+1] )*.5
                          -EB_WAVE[ii  ][jj+1][N_z-pwr_dect  ]
                          *( EB_WAVE[ii  ][jj+1][N_z-pwr_dect-1]
                            +EB_WAVE[ii  ][jj+1][N_z-pwr_dect+1] )*.5 );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x1-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[pwr_dect  ][jj+1][kk  ]
                          *( EB_WAVE[pwr_dect+1][jj+1][kk  ]
                            +EB_WAVE[pwr_dect-1][jj+1][kk  ] )*.5
                          -EB_WAVE[pwr_dect  ][jj  ][kk+1]
                          *( EB_WAVE[pwr_dect+1][jj  ][kk+1] 
                            +EB_WAVE[pwr_dect-1][jj  ][kk+1] )*.5 );
            }
        }
    } else if ( strcmp(absorber,"x2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x2-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[N_x-pwr_dect  ][jj+1][kk  ]
                          *( EB_WAVE[N_x-pwr_dect-1][jj+1][kk  ]
                            +EB_WAVE[N_x-pwr_dect+1][jj+1][kk  ] )*.5
                          -EB_WAVE[N_x-pwr_dect  ][jj  ][kk+1]
                          *( EB_WAVE[N_x-pwr_dect-1][jj  ][kk+1] 
                            +EB_WAVE[N_x-pwr_dect+1][jj  ][kk+1] )*.5 );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y1-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][pwr_dect  ][kk+1]
                          *( EB_WAVE[ii  ][pwr_dect+1][kk+1]
                            +EB_WAVE[ii  ][pwr_dect-1][kk+1] )*.5
                          -EB_WAVE[ii+1][pwr_dect  ][kk  ]
                          *( EB_WAVE[ii+1][pwr_dect+1][kk  ]
                            +EB_WAVE[ii+1][pwr_dect-1][kk  ] )*.5 );
            }
        }
    } else if ( strcmp(absorber,"y2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y2-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][N_y-pwr_dect  ][kk+1]
                          *( EB_WAVE[ii  ][N_y-pwr_dect-1][kk+1]
                            +EB_WAVE[ii  ][N_y-pwr_dect+1][kk+1] )*.5
                          -EB_WAVE[ii+1][N_y-pwr_dect  ][kk  ]
                          *( EB_WAVE[ii+1][N_y-pwr_dect-1][kk  ]
                            +EB_WAVE[ii+1][N_y-pwr_dect+1][kk  ] )*.5 );
            }
        }
    }
    
    return fabs(poynt);
} //}}}

//Print time traces to console and 
//send timetraces to grid_io for file saving
int write_timetraces(   gridConfiguration *gridCfg,
                        saveData *saveDCfg ){

    //Char values as directions to the correct folder
    char fullDir[PATH_MAX], filename_trace[PATH_MAX];

    //Create filepath to save timetraces
    sprintf(fullDir,"%s/%s/", projectPath, foldername);
    sprintf( filename_trace, "%s%s", fullDir, file_trace);

    writeConsole_timetraces( T_END, PERIOD );

    // write timetrace data into file
    writeTimetraces2ascii( T_END, PERIOD, 
                           filename_trace , timetraces );               //function in GRID_IO.C

    free2DArray( timetraces, (T_END/(int)PERIOD) );                     //free timetraces array

    return EXIT_SUCCESS;
}

//Print full timetraces to console
int writeConsole_timetraces( int T_end, double period ){

    printf( "-------------------------------------------------------------------------------------------------------------\n" );
    printf( "  T   |   poynt_z1   |   poynt_z2   |   poynt_x1   |   poynt_x2   |   poynt_y1   |   poynt_y2   |  P_out     \n" );
    printf( "------+--------------+--------------+--------------+--------------+--------------+--------------+------------\n" );
    for ( int ii=0 ; ii<(T_end/(int)period) ; ++ii )
        printf( " %4d |%13.6e |%13.6e |%13.6e |%13.6e |%13.6e |%13.6e |%13.6e\n",
                (int)timetraces[ii][1], //timetraces[ii][1],
                timetraces[ii][2], timetraces[ii][3],
                timetraces[ii][4], timetraces[ii][5],
                timetraces[ii][6], timetraces[ii][7],
                (timetraces[ii][2]+timetraces[ii][3] + timetraces[ii][4]+timetraces[ii][5] + timetraces[ii][6]+timetraces[ii][7])
              );
    printf( "-------------------------------------------------------------------------------------------------------------\n" );

    return EXIT_SUCCESS;
}

