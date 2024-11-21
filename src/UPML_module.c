#include "UPML_module.h"

static double ***DH_WAVE = NULL;
static double ***DH_WAVE_ref = NULL;

/*Initialize PML arrays functions*/
void init_UPML_fields( gridConfiguration *gridCfg ){

    DH_WAVE = allocatePMLArray(NX, NY, NZ);
    DH_WAVE_ref = allocatePMLArray(NX, NY, NZ_REF);

}

// Function to allocate memory for a 3D array with given dimensions
double ***allocatePMLArray(int N_x, int N_y, int N_z) {
    double ***array = (double ***)calloc( N_x, sizeof(double **));
    for (int i = 0; i < N_x; i++) {
        array[i] = (double **)calloc( N_y, sizeof(double *));
        for (int j = 0; j < N_y; j++) {
            array[i][j] = (double *)calloc( N_z, sizeof(double));
        }
    }
    return array;
}

/*Magnetic field UPML*/
void UPML_B_faces(  gridConfiguration *gridCfg, 
                    boundaryVariables *boundaryV,
                    double EB_WAVE[NX][NY][NZ] ){
    
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

//Boundary x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii < d_absorb ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Boundary x > Nx - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Boundary y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Boundary y > Ny - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj=NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Boundary z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Boundary z > Nz - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk= NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

}

void UPML_B_corners(gridConfiguration *gridCfg, 
                    boundaryVariables *boundaryV,
                    double EB_WAVE[NX][NY][NZ] ){

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

//Corner x, y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Corner x > Nx - d_absorb; y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Corner y > Ny - d_absorb; x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Corner x,y > N - d_absorb; z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Corner z > N - d_absorb; x,y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

//Corner x,z > N - d_absorb; y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

//Corner y,z > N - d_absorb; x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

//Corner x,y,z > N - d_absorb; 
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

}

void UPML_B_edges(  gridConfiguration *gridCfg, 
                    boundaryVariables *boundaryV,
                    double EB_WAVE[NX][NY][NZ] ){

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

//Edge x, y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Edge x < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Edge y < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Edge x, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Edge x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Edge x < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = d_absorb ; jj < NY - d_absorb - 2 ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

//Edge z < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Edge x,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

//Edge y,z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Edge y < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

//Edge z < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Edge y,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

}

void UPML_Bref_faces(   gridConfiguration *gridCfg, 
                        boundaryVariables *boundaryV,
                        double EB_WAVE[NX][NY][NZ_REF] ){
    
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

//Boundary x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii < d_absorb ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary x > Nx - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NY - d_absorb - 2 ; ii+=2) {
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary y > Ny - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary z > Nz - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

}

void UPML_Bref_corners( gridConfiguration *gridCfg, 
                        boundaryVariables *boundaryV,
                        double EB_WAVE[NX][NY][NZ_REF] ){

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

//Corner x, y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x > Nx - d_absorb; y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner y > Ny - d_absorb; x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,y > N - d_absorb; z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner z > N - d_absorb; x,y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,z > N - d_absorb; y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner y,z > N - d_absorb; x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,y,z > N - d_absorb; 
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) ) * (
                                            +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                                            -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                                            );
                EB_WAVE[ii+1][jj  ][kk+1] = Cx(ii/2)*EB_WAVE[ii+1][jj  ][kk+1] + ( 1/F2x(ii/2) )*(
                                            + F2y(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - F1y(jj/2)*dystore );

                // -dBz/DT = dEy/dx - dEx/dy
                dzstore = DH_WAVE_ref[ii+1][jj+1][kk  ];
                DH_WAVE_ref[ii+1][jj+1][kk  ] = Cx(ii/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - 1.*( 2*DT/DX/F2x(ii/2) ) * (
                                            +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                                            -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                                            );
                EB_WAVE[ii+1][jj+1][kk  ] = Cy(jj/2)*EB_WAVE[ii+1][jj+1][kk  ] + ( 1/F2y(jj/2) )*(
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

}

void UPML_Bref_edges(   gridConfiguration *gridCfg, 
                        boundaryVariables *boundaryV,
                        double EB_WAVE[NX][NY][NZ_REF] ){

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

//Edge x, y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge y < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
                
            }
        }
    }

//Edge x < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
                
            }
        }
    }

//Edge z < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge y,z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );   
            }
        }
    }

//Edge y < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );  
            }
        }
    }

//Edge z < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );    
            }
        }
    }

//Edge y,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // -dBx/DT = dEz/dy - dEy/dz
                dxstore = DH_WAVE_ref[ii  ][jj+1][kk+1];
                DH_WAVE_ref[ii  ][jj+1][kk+1] = Cy(jj/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - 1.*( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                                            -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                                            );
                EB_WAVE[ii  ][jj+1][kk+1] = Cz(kk/2)*EB_WAVE[ii  ][jj+1][kk+1] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii  ][jj+1][kk+1] - F1x(ii/2)*dxstore );

                // -dBy/DT = dEx/dz - dEz/DX
                dystore = DH_WAVE_ref[ii+1][jj  ][kk+1];
                DH_WAVE_ref[ii+1][jj  ][kk+1] = Cz(kk/2)*DH_WAVE_ref[ii+1][jj  ][kk+1] - 1.*( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii+1][jj+1][kk  ] - F1z(kk/2)*dzstore );   
            }
        }
    }

}


/*Electric field UPML*/
void UPML_E_faces(  gridConfiguration *gridCfg, 
                    boundaryVariables *boundaryV,
                    double EB_WAVE[NX][NY][NZ] ){

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

//Boundary x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii < d_absorb ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Boundary x > Nx - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = NX - d_absorb ; ii < NX - 2 ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Boundary y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Boundary y > Ny - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj= NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Boundary z < d_absorb 
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Boundary z > Nz - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk= NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

}

void UPML_E_corners(gridConfiguration *gridCfg, 
                    boundaryVariables *boundaryV,
                    double EB_WAVE[NX][NY][NZ] ){

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

//Corner x, y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Corner x > Nx - d_absorb; y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Corner y > Ny - d_absorb; x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Corner x,y > N - d_absorb; z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Corner z > N - d_absorb; x,y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

//Corner x,z > N - d_absorb; y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

//Corner y,z > N - d_absorb; x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

//Corner x,y,z > N - d_absorb; 
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

}

void UPML_E_edges(  gridConfiguration *gridCfg, 
                    boundaryVariables *boundaryV,
                    double EB_WAVE[NX][NY][NZ] ){

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

//Corner x, y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Edge x < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Edge y < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Edge x, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ - d_absorb - 2 ; kk+=2) {
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

//Edge x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Edge x < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

//Edge z < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Edge x,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

//Edge y,z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Edge y < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

//Edge z < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
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

//Edge y,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = NZ - d_absorb ; kk < NZ-2 ; kk+=2) {
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

}

void UPML_Eref_faces(   gridConfiguration *gridCfg, 
                        boundaryVariables *boundaryV,
                        double EB_WAVE[NX][NY][NZ_REF] ){

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

//Boundary x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii < d_absorb ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary x > Nx - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = NX - d_absorb ; ii < NX - 2 ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary y > Ny - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Boundary z > Nz - d_absorb - 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

}

void UPML_Eref_corners( gridConfiguration *gridCfg, 
                        boundaryVariables *boundaryV,
                        double EB_WAVE[NX][NY][NZ_REF]){

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

//Corner x, y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x > Nx - d_absorb; y, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner y > Ny - d_absorb; x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,y > N - d_absorb; z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner z > N - d_absorb; x,y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,z > N - d_absorb; y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner y,z > N - d_absorb; x < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Corner x,y,z > N - d_absorb; 
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

}

void UPML_Eref_edges(   gridConfiguration *gridCfg, 
                        boundaryVariables *boundaryV,
                        double EB_WAVE[NX][NY][NZ_REF] ){

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

//Edge x, y < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge y < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = d_absorb ; kk <= NZ_REF - d_absorb - 2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x, z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
                
            }
        }
    }

//Edge x < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii < d_absorb ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
                
            }
        }
    }

//Edge z < d_absorb + 2, x > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge x,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = NX - d_absorb ; ii < NX-2 ; ii+=2) {                
        for (jj = d_absorb ; jj <= NY - d_absorb - 2 ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );
            }
        }
    }

//Edge y,z < d_absorb + 2
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );  
            }
        }
    }

//Edge y < d_absorb + 2, z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj=2 ; jj < d_absorb ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );  
            }
        }
    }

//Edge z < d_absorb + 2, y > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii= d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk=2 ; kk < d_absorb ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );   
            }
        }
    }

//Edge y,z > N - d_absorb
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii = d_absorb ; ii <= NX - d_absorb - 2 ; ii+=2) {                
        for (jj = NY - d_absorb ; jj < NY-2 ; jj+=2) {
            for (kk = NZ_REF - d_absorb ; kk < NZ_REF-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                dxstore = DH_WAVE_ref[ii+1][jj  ][kk  ];
                DH_WAVE_ref[ii+1][jj  ][kk  ] = Cy(jj/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] + ( 2*DT/DX/F2y(jj/2) )*(
                                            +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                                            -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                                            );
                EB_WAVE[ii+1][jj  ][kk  ] = Cz(kk/2)*EB_WAVE[ii+1][jj  ][kk  ] + ( 1/F2z(kk/2) )*(
                                            + F2x(ii/2)*DH_WAVE_ref[ii+1][jj  ][kk  ] - F1x(ii/2)*dxstore );

                // dEy/dt = (dBx/dz - dBz/dx)
                dystore = DH_WAVE_ref[ii  ][jj+1][kk  ];
                DH_WAVE_ref[ii  ][jj+1][kk  ] = Cz(kk/2)*DH_WAVE_ref[ii  ][jj+1][kk  ] + ( 2*DT/DX/F2z(kk/2) )*(
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
                                            + F2z(kk/2)*DH_WAVE_ref[ii  ][jj  ][kk+1] - F1z(kk/2)*dzstore );    
            }
        }
    }

}

