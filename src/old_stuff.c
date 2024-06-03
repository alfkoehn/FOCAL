/* This file is meant to be used as a repository for old and no longer used
 * functions. In an ideal world, there would be no need to ever look at these
 * functions again. We are not living in an ideal world though, therefore I
 * decided to move all old functions in this file. Note that some of these 
 * functions were considered as wrong when they were moved into this file.
 */


// taken from src/power_calc.c on 2024-03-06
double calc_power_EE_1( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                       int pwr_dect, char absorber[],
                       double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] ) {
//{{{
    // sum-up the power (squared E-field) in the absorbers

    size_t
        ii, jj, kk;
    double
        power;

    power   = .0;

    // Bx: even-odd-odd
    // By: odd-even-odd
    // Bz: odd-odd-even
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    if ( strcmp(absorber,"ref_z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:power)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                power += ( pow(EB_WAVE_ref[ii+1][jj  ][pwr_dect  ], 2)
                          +pow(EB_WAVE_ref[ii  ][jj+1][pwr_dect  ], 2)
                          +pow(EB_WAVE_ref[ii  ][jj  ][pwr_dect+1], 2) );
            }
        }
    } else if ( strcmp(absorber,"z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:power)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                power += ( pow((EB_WAVE[ii+1][jj  ][pwr_dect  ] - EB_WAVE_ref[ii+1][jj  ][pwr_dect  ]), 2)
                          +pow((EB_WAVE[ii  ][jj+1][pwr_dect  ] - EB_WAVE_ref[ii  ][jj+1][pwr_dect  ]), 2)
                          +pow((EB_WAVE[ii  ][jj  ][pwr_dect+1] - EB_WAVE_ref[ii  ][jj  ][pwr_dect+1]), 2) );
            }
        }
    } else if ( strcmp(absorber,"z2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:power)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z2-plane
                power += ( pow(EB_WAVE[ii+1][jj  ][N_z-pwr_dect  ], 2)
                          +pow(EB_WAVE[ii  ][jj+1][N_z-pwr_dect  ], 2)
                          +pow(EB_WAVE[ii  ][jj  ][N_z-pwr_dect-1], 2) );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:power)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x1-plane
                power += ( pow(EB_WAVE[pwr_dect+1][jj  ][kk  ], 2)
                          +pow(EB_WAVE[pwr_dect  ][jj+1][kk  ], 2)
                          +pow(EB_WAVE[pwr_dect  ][jj  ][kk+1], 2) );
            }
        }
    } else if ( strcmp(absorber,"x2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:power)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x2-plane
                power += ( pow(EB_WAVE[N_x-pwr_dect-1][jj  ][kk  ], 2)
                          +pow(EB_WAVE[N_x-pwr_dect  ][jj+1][kk  ], 2)
                          +pow(EB_WAVE[N_x-pwr_dect  ][jj  ][kk+1], 2) );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:power)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y1-plane
                power += ( pow(EB_WAVE[ii+1][pwr_dect  ][kk  ], 2)
                          +pow(EB_WAVE[ii  ][pwr_dect+1][kk  ], 2)
                          +pow(EB_WAVE[ii  ][pwr_dect  ][kk+1], 2) );
            }
        }
    } else if ( strcmp(absorber,"y2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:power)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y2-plane
                power += ( pow(EB_WAVE[ii+1][N_y-pwr_dect  ][kk  ], 2)
                          +pow(EB_WAVE[ii  ][N_y-pwr_dect+1][kk  ], 2)
                          +pow(EB_WAVE[ii  ][N_y-pwr_dect  ][kk+1], 2) );
            }
        }
    }
    
    return fabs(power);
} //}}}

// taken from include/power_calc.h on 2024-03-06
double calc_power_EE_1( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                       int pwr_dect, char absorber[],
                       double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] );


