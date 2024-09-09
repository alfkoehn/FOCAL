#include "focal.h"
#include "focal-struct.h"

int advance_J( gridConfiguration *gridCfg, 
               double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
               double J_B0[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz],
               double n_e[gridCfg->Nx/2][gridCfg->Ny/2][gridCfg->Nz/2] ) { 
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
    // B_new = B_old - nabla x E

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
    // E_new = E_old + c^2*nablaxB - 1/epsilon_0*J

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


