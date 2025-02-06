#include "UPML_module.h"

static double ***DH_WAVE = NULL;
static double ***DH_WAVE_ref = NULL;


/*Initialize PML arrays functions*/
void init_UPML_fields( gridConfiguration *gridCfg ){
    //{{{

    DH_WAVE = allocate3DArray(NX, NY, NZ);
    DH_WAVE_ref = allocate3DArray(NX, NY, NZ_REF);

}//}}}


int free_UPML_memory( gridConfiguration *gridCfg ){
    //{{{

    free3DArray(DH_WAVE, NX, NY);
    free3DArray(DH_WAVE_ref, NX, NY);

    return EXIT_SUCCESS;
}//}}}


/*Magnetic field UPML*/
void UPML_B_faces(  gridConfiguration *gridCfg, 
                    boundaryVariables *boundaryV,
                    double EB_WAVE[NX][NY][NZ] ){
    //{{{
    
    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Boundary x < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary x > Nx - D_ABSORB - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary y > Ny - D_ABSORB - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj=NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary z > Nz - D_ABSORB - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk= NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

}//}}}


void UPML_B_corners(gridConfiguration *gridCfg, 
                    boundaryVariables *boundaryV,
                    double EB_WAVE[NX][NY][NZ] ){
    //{{{

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Corner x, y, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x > Nx - D_ABSORB; y, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner y > Ny - D_ABSORB; x, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,y > N - D_ABSORB; z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner z > N - D_ABSORB; x,y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,z > N - D_ABSORB; y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner y,z > N - D_ABSORB; x < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,y,z > N - D_ABSORB; 
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

}//}}}


void UPML_B_edges(  gridConfiguration *gridCfg, 
                    boundaryVariables *boundaryV,
                    double EB_WAVE[NX][NY][NZ] ){
    //{{{

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Edge x, y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x < D_ABSORB + 2, y > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge y < D_ABSORB + 2, x > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x, y > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
                
            }
        }
    }

//Edge x < D_ABSORB + 2, z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = D_ABSORB ; jj < NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
                
            }
        }
    }

//Edge z < D_ABSORB + 2, x > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x,z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge y,z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );   
            }
        }
    }

//Edge y < D_ABSORB + 2, z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );   
            }
        }
    }

//Edge z < D_ABSORB + 2, y > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );    
            }
        }
    }

//Edge y,z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE[ii  ][jj+1][kk+1];
                DH_WAVE[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE[ii+1][jj  ][kk+1];
                DH_WAVE[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE[ii+1][jj+1][kk  ];
                DH_WAVE[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );   
            }
        }
    }

}//}}}


void UPML_Bref_faces(   gridConfiguration *gridCfg, 
                        boundaryVariables *boundaryV,
                        double EB_WAVE[NX][NY][NZ_REF] ){
    //{{{
    
    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Boundary x < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary x > Nx - D_ABSORB - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NY - D_ABSORB - 2 ; ii+=2) {
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary y > Ny - D_ABSORB - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary z > Nz - D_ABSORB - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

}//}}}


void UPML_Bref_corners( gridConfiguration *gridCfg, 
                        boundaryVariables *boundaryV,
                        double EB_WAVE[NX][NY][NZ_REF] ){
    //{{{

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Corner x, y, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x > Nx - D_ABSORB; y, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner y > Ny - D_ABSORB; x, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x,y > N - D_ABSORB; z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner z > N - D_ABSORB; x,y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x,z > N - D_ABSORB; y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner y,z > N - D_ABSORB; x < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x,y,z > N - D_ABSORB; 
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

}//}}}


void UPML_Bref_edges(   gridConfiguration *gridCfg, 
                        boundaryVariables *boundaryV,
                        double EB_WAVE[NX][NY][NZ_REF] ){
    //{{{

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Edge x, y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x < D_ABSORB + 2, y > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge y < D_ABSORB + 2, x > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x, y > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
                
            }
        }
    }

//Edge x < D_ABSORB + 2, z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
                
            }
        }
    }

//Edge z < D_ABSORB + 2, x > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x,z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge y,z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );   
            }
        }
    }

//Edge y < D_ABSORB + 2, z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );  
            }
        }
    }

//Edge z < D_ABSORB + 2, y > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );   
            }
        }
    }

//Edge y,z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Czr(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Czr(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1zr(kk/2)*dzstore );   
            }
        }
    }

}//}}}


/*Electric field UPML*/
void UPML_E_faces(  gridConfiguration *gridCfg, 
                    boundaryVariables *boundaryV,
                    double EB_WAVE[NX][NY][NZ] ){
    //{{{

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Boundary x < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary x > Nx - D_ABSORB - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = NX - D_ABSORB ; ii < NX - 2 ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary y > Ny - D_ABSORB - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj= NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary z < D_ABSORB 
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary z > Nz - D_ABSORB - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk= NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

}//}}}


void UPML_E_corners(gridConfiguration *gridCfg, 
                    boundaryVariables *boundaryV,
                    double EB_WAVE[NX][NY][NZ] ){
    //{{{

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Corner x, y, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x > Nx - D_ABSORB; y, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner y > Ny - D_ABSORB; x, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,y > N - D_ABSORB; z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner z > N - D_ABSORB; x,y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,z > N - D_ABSORB; y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner y,z > N - D_ABSORB; x < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,y,z > N - D_ABSORB; 
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

}//}}}


void UPML_E_edges(  gridConfiguration *gridCfg, 
                    boundaryVariables *boundaryV,
                    double EB_WAVE[NX][NY][NZ] ){
    //{{{

    // DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Corner x, y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x < D_ABSORB + 2, y > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge y < D_ABSORB + 2, x > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x, y > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x < D_ABSORB + 2, z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge z < D_ABSORB + 2, x > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x,z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge y,z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge y < D_ABSORB + 2, z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge z < D_ABSORB + 2, y > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );   
            }
        }
    }

//Edge y,z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = NZ - D_ABSORB ; kk < NZ-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE[ii+1][jj  ][kk  ];
                DH_WAVE[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) ) * (
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE[ii  ][jj+1][kk  ];
                DH_WAVE[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE[ii  ][jj  ][kk+1];
                DH_WAVE[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

}//}}}


void UPML_Eref_faces(   gridConfiguration *gridCfg, 
                        boundaryVariables *boundaryV,
                        double EB_WAVE[NX][NY][NZ_REF] ){
    //{{{

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Boundary x < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary x > Nx - D_ABSORB - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = NX - D_ABSORB ; ii < NX - 2 ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary y > Ny - D_ABSORB - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Boundary z > Nz - D_ABSORB - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

}//}}}


void UPML_Eref_corners( gridConfiguration *gridCfg, 
                        boundaryVariables *boundaryV,
                        double EB_WAVE[NX][NY][NZ_REF]){
    //{{{

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Corner x, y, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x > Nx - D_ABSORB; y, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner y > Ny - D_ABSORB; x, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x,y > N - D_ABSORB; z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner z > N - D_ABSORB; x,y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x,z > N - D_ABSORB; y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner y,z > N - D_ABSORB; x < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Corner x,y,z > N - D_ABSORB; 
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

}//}}}


void UPML_Eref_edges(   gridConfiguration *gridCfg, 
                        boundaryVariables *boundaryV,
                        double EB_WAVE[NX][NY][NZ_REF] ){
    //{{{

    //DH_WAVE:                  EB_WAVE:
    // Dx: odd-even-even        Ex:odd-even-even
    // Dy: even-odd-even        Ey: even-odd-even
    // Dz: even-even-odd        Ez:even-even-odd
    // Hx: even-odd-odd         Bx: even-odd-odd
    // Hy: odd-even-odd         By: odd-even-odd
    // Hz: odd-odd-even         Bz: odd-odd-even        

    size_t
        ii, jj, kk;
    double
        dxstore, dystore, dzstore;

//Edge x, y < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x < D_ABSORB + 2, y > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge y < D_ABSORB + 2, x > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x, y > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = D_ABSORB ; kk <= NZ_REF - D_ABSORB - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x, z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
                
            }
        }
    }

//Edge x < D_ABSORB + 2, z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < D_ABSORB ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
                
            }
        }
    }

//Edge z < D_ABSORB + 2, x > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge x,z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - D_ABSORB ; ii < NX-2 ; ii+=2) {                
        for (jj = D_ABSORB ; jj <= NY - D_ABSORB - 2 ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );
            }
        }
    }

//Edge y,z < D_ABSORB + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );  
            }
        }
    }

//Edge y < D_ABSORB + 2, z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj=2 ; jj < D_ABSORB ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );  
            }
        }
    }

//Edge z < D_ABSORB + 2, y > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < D_ABSORB ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );   
            }
        }
    }

//Edge y,z > N - D_ABSORB
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = D_ABSORB ; ii <= NX - D_ABSORB - 2 ; ii+=2) {                
        for (jj = NY - D_ABSORB ; jj < NY-2 ; jj+=2) {
            for (kk = NZ_REF - D_ABSORB ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Czr(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2zr(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Czr(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2zr(kk/2) )*(
                                            +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                                            -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk  ] = Cx(ii/2)*EB_WAVE[ii  ][jj+1][kk  ] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] - F1y(jj/2)*dystore );

                // dEz/dt = (dBy/dx - dBx/dy)
                dzstore = DH_WAVE_ref[ii  ][jj  ][kk+1];
                DH_WAVE_ref[ii  ][jj  ][kk+1] = Cx(ii/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] + ( 2*DT/DX/F2x(ii/2) )*(
                                            +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                                            );
                EB_WAVE[ii  ][jj  ][kk+1] = Cy(jj/2)*EB_WAVE[ii  ][jj  ][kk+1] + ( 1/F2y(jj/2) )*(
                                            + F2zr(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1zr(kk/2)*dzstore );    
            }
        }
    }

}//}}}

