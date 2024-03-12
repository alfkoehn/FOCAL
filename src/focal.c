#include <math.h>
#include <stdlib.h>

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


