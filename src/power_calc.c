#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "focal-struct.h"
#include "focal.h"
#include "power_calc.h"

double calc_poynt_4( gridConfiguration *gridCfg, 
                     int pwr_dect, char absorber[],
                     double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                     double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref] ) {
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
        for (ii=pwr_dect ; ii<=(gridCfg->Nx-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<=(gridCfg->Ny-pwr_dect-2) ; jj+=2) {
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
        for (ii=pwr_dect ; ii<=(gridCfg->Nx-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<=(gridCfg->Ny-pwr_dect-2) ; jj+=2) {
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
        for (ii=pwr_dect ; ii<=(gridCfg->Nx-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<=(gridCfg->Ny-pwr_dect-2) ; jj+=2) {
                // z2-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE[ii+1][jj  ][gridCfg->Nz-pwr_dect  ]
                          *EB_WAVE[ii+1][jj  ][gridCfg->Nz-pwr_dect-1]
                          -EB_WAVE[ii  ][jj+1][gridCfg->Nz-pwr_dect  ]
                          *EB_WAVE[ii  ][jj+1][gridCfg->Nz-pwr_dect-1] );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(gridCfg->Ny-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(gridCfg->Nz-pwr_dect-2) ; kk+=2) {
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
        for (jj=pwr_dect ; jj<=(gridCfg->Ny-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(gridCfg->Nz-pwr_dect-2) ; kk+=2) {
                // x2-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[gridCfg->Nx-pwr_dect  ][jj+1][kk  ]
                          *EB_WAVE[gridCfg->Nx-pwr_dect-1][jj+1][kk  ]
                          -EB_WAVE[gridCfg->Nx-pwr_dect  ][jj  ][kk+1]
                          *EB_WAVE[gridCfg->Nx-pwr_dect-1][jj  ][kk+1] );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(gridCfg->Nx-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(gridCfg->Nz-pwr_dect-2) ; kk+=2) {
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
        for (ii=pwr_dect ; ii<=(gridCfg->Nx-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(gridCfg->Nz-pwr_dect-2) ; kk+=2) {
                // y2-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][gridCfg->Ny-pwr_dect  ][kk+1]
                          *EB_WAVE[ii  ][gridCfg->Ny-pwr_dect-1][kk+1]
                          -EB_WAVE[ii+1][gridCfg->Ny-pwr_dect  ][kk  ]
                          *EB_WAVE[ii+1][gridCfg->Ny-pwr_dect-1][kk  ] );
            }
        }
    }
    
    return fabs(poynt);
} //}}}


double calc_poynt_5( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     int pwr_dect, char absorber[],
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
                     int pwr_dect, char absorber[],
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
                     int pwr_dect, char absorber[],
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


