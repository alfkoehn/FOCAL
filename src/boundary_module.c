#include "boundary_module.h"

static double ***E_Xdir_OLD = NULL;
static double ***E_Ydir_OLD = NULL;
static double ***E_Zdir_OLD = NULL;
static double ***E_Xdir_OLD_ref = NULL;
static double ***E_Ydir_OLD_ref = NULL;
static double ***E_Zdir_OLD_ref = NULL;

void init_boundary(gridConfiguration *gridCfg, boundaryVariables *boundaryV){

    if(sel_boundary == 1){

        eco = 10./(double)(period);

    }
    else if (sel_boundary == 2){
        
        E_Xdir_OLD = allocateBoundaryArray(8, NY, NZ);
        E_Ydir_OLD = allocateBoundaryArray(NX, 8, NZ);
        E_Zdir_OLD = allocateBoundaryArray(NX, NY, 8);

        E_Xdir_OLD_ref = allocateBoundaryArray(8, NY, NZ_REF);
        E_Ydir_OLD_ref = allocateBoundaryArray(NX, 8, NZ_REF);
        E_Zdir_OLD_ref = allocateBoundaryArray(NX, NY, 8);

    }
    else if(sel_boundary == 3){
        
        /*Initialize UPML parameters*/
        ALLOC_1D(boundaryV->F1x, NX/2, double);
        ALLOC_1D(boundaryV->F1y, NY/2, double);
        ALLOC_1D(boundaryV->F1z, NZ/2, double);
        ALLOC_1D(boundaryV->F2x, NX/2, double);
        ALLOC_1D(boundaryV->F2y, NY/2, double);
        ALLOC_1D(boundaryV->F2z, NZ/2, double);
        ALLOC_1D(boundaryV->Cx, NX/2, double);
        ALLOC_1D(boundaryV->Cy, NY/2, double);
        ALLOC_1D(boundaryV->Cz, NZ/2, double);

        /*UPML ref parameters*/
        ALLOC_1D(boundaryV->F1zr, NZ_REF/2, double);
        ALLOC_1D(boundaryV->F2zr, NZ_REF/2, double);
        ALLOC_1D(boundaryV->Czr, NZ_REF/2, double);

        init_UPML_parameters( gridCfg, boundaryV );
        init_UPML_fields( gridCfg );

    }
   
}

// Function to allocate memory for a 3D array with given dimensions
double ***allocateBoundaryArray(int N_x, int N_y, int N_z) {
    double ***array = (double ***)calloc( N_x, sizeof(double **));
    for (int i = 0; i < N_x; i++) {
        array[i] = (double **)calloc( N_y, sizeof(double *));
        for (int j = 0; j < N_y; j++) {
            array[i][j] = (double *)calloc( N_z, sizeof(double));
        }
    }
    return array;
}

/*Apply boundary on time evolution*/
void advance_boundary(  gridConfiguration *gridCfg, boundaryVariables *boundaryV, 
                        double EB_WAVE[NX][NY][NZ], double EB_WAVE_ref[NX][NY][NZ_REF]){

    if(sel_boundary == 1){

        apply_absorber( gridCfg, boundaryV, EB_WAVE);
        apply_absorber_ref(gridCfg, boundaryV, EB_WAVE_ref);

    }else if (sel_boundary == 2){

        abc_Mur_1st( gridCfg, "x1x2y1y2z1z2", EB_WAVE );
        abc_Mur_1st_ref( gridCfg, EB_WAVE_ref );
                         
        abc_Mur_saveOldE_xdir(    gridCfg, EB_WAVE );
        abc_Mur_saveOldE_ydir(    gridCfg, EB_WAVE );
        abc_Mur_saveOldE_zdir(    gridCfg, EB_WAVE );
        abc_Mur_saveOldEref_xdir( gridCfg, EB_WAVE_ref );
        abc_Mur_saveOldEref_ydir( gridCfg, EB_WAVE_ref );
        abc_Mur_saveOldEref_zdir( gridCfg, EB_WAVE_ref );
        
    }else if (sel_boundary == 3){

        UPML_B_faces(   gridCfg, boundaryV, EB_WAVE );
        UPML_B_corners( gridCfg, boundaryV, EB_WAVE );
        UPML_B_edges(   gridCfg, boundaryV, EB_WAVE );

        UPML_Bref_faces(    gridCfg, boundaryV, EB_WAVE_ref );
        UPML_Bref_corners(  gridCfg, boundaryV, EB_WAVE_ref );
        UPML_Bref_edges(    gridCfg, boundaryV, EB_WAVE_ref );

        UPML_E_faces(   gridCfg, boundaryV, EB_WAVE );
        UPML_E_corners( gridCfg, boundaryV, EB_WAVE );
        UPML_E_edges(   gridCfg, boundaryV, EB_WAVE );

        UPML_Eref_faces(    gridCfg, boundaryV, EB_WAVE_ref );
        UPML_Eref_corners(  gridCfg, boundaryV, EB_WAVE_ref );
        UPML_Eref_edges(    gridCfg, boundaryV, EB_WAVE_ref );

    }

}


/*Boundary Functions*/
/*ABC functions*/
int apply_absorber( gridConfiguration *gridCfg, 
                    boundaryVariables *boundaryV, 
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
                        boundaryVariables *boundaryV, 
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


int apply_absorber_v2(  size_t N_x, size_t N_y, size_t N_z, int D_absorb, 
                        boundaryVariables *boundaryV, 
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

/*Mur boundary functions*/
int abc_Mur_saveOldE_xdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[NX][NY][NZ] ) {
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
            E_Xdir_OLD[0+1][jj  ][kk  ]  = EB_WAVE[0+offset+1][jj  ][kk  ];
            E_Xdir_OLD[2+1][jj  ][kk  ]  = EB_WAVE[2+offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_Xdir_OLD[0  ][jj+1][kk  ]  = EB_WAVE[0+offset  ][jj+1][kk  ];
            E_Xdir_OLD[2  ][jj+1][kk  ]  = EB_WAVE[2+offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_Xdir_OLD[0  ][jj  ][kk+1]  = EB_WAVE[0+offset  ][jj  ][kk+1];
            E_Xdir_OLD[2  ][jj  ][kk+1]  = EB_WAVE[2+offset  ][jj  ][kk+1];

            // store values at x=NX-1 and x=NX-2
            // Ex: odd-even-even
            E_Xdir_OLD[4+1][jj  ][kk  ]  = EB_WAVE[NX-4-offset+1][jj  ][kk  ];
            E_Xdir_OLD[6+1][jj  ][kk  ]  = EB_WAVE[NX-2-offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_Xdir_OLD[4  ][jj+1][kk  ]  = EB_WAVE[NX-4-offset  ][jj+1][kk  ];
            E_Xdir_OLD[6  ][jj+1][kk  ]  = EB_WAVE[NX-2-offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_Xdir_OLD[4  ][jj  ][kk+1]  = EB_WAVE[NX-4-offset  ][jj  ][kk+1];
            E_Xdir_OLD[6  ][jj  ][kk+1]  = EB_WAVE[NX-2-offset  ][jj  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}

int abc_Mur_saveOldE_ydir( gridConfiguration *gridCfg, 
                           double EB_WAVE[NX][NY][NZ] ) {
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
            E_Ydir_OLD[ii+1][0  ][kk  ]  = EB_WAVE[ii+1][0+offset  ][kk  ];
            E_Ydir_OLD[ii+1][2  ][kk  ]  = EB_WAVE[ii+1][2+offset  ][kk  ];
            // Ey: even-odd-even
            E_Ydir_OLD[ii  ][0+1][kk  ]  = EB_WAVE[ii  ][0+offset+1][kk  ];
            E_Ydir_OLD[ii  ][2+1][kk  ]  = EB_WAVE[ii  ][2+offset+1][kk  ];
            // Ez: even-even-odd
            E_Ydir_OLD[ii  ][0  ][kk+1]  = EB_WAVE[ii  ][0+offset  ][kk+1];
            E_Ydir_OLD[ii  ][2  ][kk+1]  = EB_WAVE[ii  ][2+offset  ][kk+1];

            // store values at x=NX-1 and x=NX-2
            // Ex: odd-even-even
            E_Ydir_OLD[ii+1][4  ][kk  ]  = EB_WAVE[ii+1][NY-4-offset  ][kk  ];
            E_Ydir_OLD[ii+1][6  ][kk  ]  = EB_WAVE[ii+1][NY-2-offset  ][kk  ];
            // Ey: even-odd-even
            E_Ydir_OLD[ii  ][4+1][kk  ]  = EB_WAVE[ii  ][NY-4-offset+1][kk  ];
            E_Ydir_OLD[ii  ][6+1][kk  ]  = EB_WAVE[ii  ][NY-2-offset+1][kk  ];
            // Ez: even-even-odd
            E_Ydir_OLD[ii  ][4  ][kk+1]  = EB_WAVE[ii  ][NY-4-offset  ][kk+1];
            E_Ydir_OLD[ii  ][6  ][kk+1]  = EB_WAVE[ii  ][NY-2-offset  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}

int abc_Mur_saveOldE_zdir(  gridConfiguration *gridCfg, 
                            double EB_WAVE[NX][NY][NZ] ) {
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
            E_Zdir_OLD[ii+1][jj  ][0  ]  = EB_WAVE[ii+1][jj  ][0+offset  ];
            E_Zdir_OLD[ii+1][jj  ][2  ]  = EB_WAVE[ii+1][jj  ][2+offset  ];
            // Ey: even-odd-even
            E_Zdir_OLD[ii  ][jj+1][0  ]  = EB_WAVE[ii  ][jj+1][0+offset  ];
            E_Zdir_OLD[ii  ][jj+1][2  ]  = EB_WAVE[ii  ][jj+1][2+offset  ];
            // Ez: even-even-odd
            E_Zdir_OLD[ii  ][jj  ][0+1]  = EB_WAVE[ii  ][jj  ][0+offset+1];
            E_Zdir_OLD[ii  ][jj  ][2+1]  = EB_WAVE[ii  ][jj  ][2+offset+1];

            // store values at z=NZ-1 and z=NZ-2
            // Ex: odd-even-even
            E_Zdir_OLD[ii+1][jj  ][4  ]  = EB_WAVE[ii+1][jj  ][NZ-4-offset  ];
            E_Zdir_OLD[ii+1][jj  ][6  ]  = EB_WAVE[ii+1][jj  ][NZ-2-offset  ];
            // Ey: even-odd-even
            E_Zdir_OLD[ii  ][jj+1][4  ]  = EB_WAVE[ii  ][jj+1][NZ-4-offset  ];
            E_Zdir_OLD[ii  ][jj+1][6  ]  = EB_WAVE[ii  ][jj+1][NZ-2-offset  ];
            // Ez: even-even-odd
            E_Zdir_OLD[ii  ][jj  ][4+1]  = EB_WAVE[ii  ][jj  ][NZ-4-offset+1];
            E_Zdir_OLD[ii  ][jj  ][6+1]  = EB_WAVE[ii  ][jj  ][NZ-2-offset+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}

int abc_Mur_saveOldEref_xdir(  gridConfiguration *gridCfg, 
                                double EB_WAVE_ref[NX][NY][NZ_REF] ) {
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
            E_Xdir_OLD_ref[0+1][jj  ][kk  ]  = EB_WAVE_ref[0+offset+1][jj  ][kk  ];
            E_Xdir_OLD_ref[2+1][jj  ][kk  ]  = EB_WAVE_ref[2+offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_Xdir_OLD_ref[0  ][jj+1][kk  ]  = EB_WAVE_ref[0+offset  ][jj+1][kk  ];
            E_Xdir_OLD_ref[2  ][jj+1][kk  ]  = EB_WAVE_ref[2+offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_Xdir_OLD_ref[0  ][jj  ][kk+1]  = EB_WAVE_ref[0+offset  ][jj  ][kk+1];
            E_Xdir_OLD_ref[2  ][jj  ][kk+1]  = EB_WAVE_ref[2+offset  ][jj  ][kk+1];

            // store values at x=NX-1 and x=NX-2
            // Ex: odd-even-even
            E_Xdir_OLD_ref[4+1][jj  ][kk  ]  = EB_WAVE_ref[NX-4-offset+1][jj  ][kk  ];
            E_Xdir_OLD_ref[6+1][jj  ][kk  ]  = EB_WAVE_ref[NX-2-offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_Xdir_OLD_ref[4  ][jj+1][kk  ]  = EB_WAVE_ref[NX-4-offset  ][jj+1][kk  ];
            E_Xdir_OLD_ref[6  ][jj+1][kk  ]  = EB_WAVE_ref[NX-2-offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_Xdir_OLD_ref[4  ][jj  ][kk+1]  = EB_WAVE_ref[NX-4-offset  ][jj  ][kk+1];
            E_Xdir_OLD_ref[6  ][jj  ][kk+1]  = EB_WAVE_ref[NX-2-offset  ][jj  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldEref_ydir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[NX][NY][NZ_REF] ) {
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
            E_Ydir_OLD_ref[ii+1][0  ][kk  ]  = EB_WAVE_ref[ii+1][0+offset  ][kk  ];
            E_Ydir_OLD_ref[ii+1][2  ][kk  ]  = EB_WAVE_ref[ii+1][2+offset  ][kk  ];
            // Ey: even-odd-even
            E_Ydir_OLD_ref[ii  ][0+1][kk  ]  = EB_WAVE_ref[ii  ][0+offset+1][kk  ];
            E_Ydir_OLD_ref[ii  ][2+1][kk  ]  = EB_WAVE_ref[ii  ][2+offset+1][kk  ];
            // Ez: even-even-odd
            E_Ydir_OLD_ref[ii  ][0  ][kk+1]  = EB_WAVE_ref[ii  ][0+offset  ][kk+1];
            E_Ydir_OLD_ref[ii  ][2  ][kk+1]  = EB_WAVE_ref[ii  ][2+offset  ][kk+1];

            // store values at x=NX-1 and x=NX-2
            // Ex: odd-even-even
            E_Ydir_OLD_ref[ii+1][4  ][kk  ]  = EB_WAVE_ref[ii+1][NY-4-offset  ][kk  ];
            E_Ydir_OLD_ref[ii+1][6  ][kk  ]  = EB_WAVE_ref[ii+1][NY-2-offset  ][kk  ];
            // Ey: even-odd-even
            E_Ydir_OLD_ref[ii  ][4+1][kk  ]  = EB_WAVE_ref[ii  ][NY-4-offset+1][kk  ];
            E_Ydir_OLD_ref[ii  ][6+1][kk  ]  = EB_WAVE_ref[ii  ][NY-2-offset+1][kk  ];
            // Ez: even-even-odd
            E_Ydir_OLD_ref[ii  ][4  ][kk+1]  = EB_WAVE_ref[ii  ][NY-4-offset  ][kk+1];
            E_Ydir_OLD_ref[ii  ][6  ][kk+1]  = EB_WAVE_ref[ii  ][NY-2-offset  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldEref_zdir(   gridConfiguration *gridCfg, 
                                double EB_WAVE_ref[NX][NY][NZ_REF] ) {
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
            E_Zdir_OLD_ref[ii+1][jj  ][0  ]  = EB_WAVE_ref[ii+1][jj  ][0+offset  ];
            E_Zdir_OLD_ref[ii+1][jj  ][2  ]  = EB_WAVE_ref[ii+1][jj  ][2+offset  ];
            // Ey: even-odd-even
            E_Zdir_OLD_ref[ii  ][jj+1][0  ]  = EB_WAVE_ref[ii  ][jj+1][0+offset  ];
            E_Zdir_OLD_ref[ii  ][jj+1][2  ]  = EB_WAVE_ref[ii  ][jj+1][2+offset  ];
            // Ez: even-even-odd
            E_Zdir_OLD_ref[ii  ][jj  ][0+1]  = EB_WAVE_ref[ii  ][jj  ][0+offset+1];
            E_Zdir_OLD_ref[ii  ][jj  ][2+1]  = EB_WAVE_ref[ii  ][jj  ][2+offset+1];

            // store values at z=NZ-1 and z=NZ-2
            // Ex: odd-even-even
            E_Zdir_OLD_ref[ii+1][jj  ][4  ]  = EB_WAVE_ref[ii+1][jj  ][NZ_REF-4-offset  ];
            E_Zdir_OLD_ref[ii+1][jj  ][6  ]  = EB_WAVE_ref[ii+1][jj  ][NZ_REF-2-offset  ];
            // Ey: even-odd-even
            E_Zdir_OLD_ref[ii  ][jj+1][4  ]  = EB_WAVE_ref[ii  ][jj+1][NZ_REF-4-offset  ];
            E_Zdir_OLD_ref[ii  ][jj+1][6  ]  = EB_WAVE_ref[ii  ][jj+1][NZ_REF-2-offset  ];
            // Ez: even-even-odd
            E_Zdir_OLD_ref[ii  ][jj  ][4+1]  = EB_WAVE_ref[ii  ][jj  ][NZ_REF-4-offset+1];
            E_Zdir_OLD_ref[ii  ][jj  ][6+1]  = EB_WAVE_ref[ii  ][jj  ][NZ_REF-2-offset+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_1st( gridConfiguration *gridCfg, 
                 char absorber[],
                 double EB_WAVE[NX][NY][NZ] ) {
//{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk,
        offset;             // refers to EB_WAVE only

    double
        cnst;

    cnst    = (DT-DX)/(DT+DX);
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
                EB_WAVE[offset+0+1][jj  ][kk  ] = E_Xdir_OLD[2+1][jj  ][kk  ]
                    + cnst * (    EB_WAVE[offset+2+1][jj  ][kk  ] 
                              -E_Xdir_OLD[0+1       ][jj  ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[offset+0  ][jj+1][kk  ] = E_Xdir_OLD[2  ][jj+1][kk  ]
                    + cnst * (    EB_WAVE[offset+2  ][jj+1][kk  ] 
                              -E_Xdir_OLD[0         ][jj+1][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[offset+0  ][jj  ][kk+1] = E_Xdir_OLD[2  ][jj  ][kk+1]
                    + cnst * (    EB_WAVE[offset+2  ][jj  ][kk+1] 
                              -E_Xdir_OLD[0         ][jj  ][kk+1] );
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
                EB_WAVE[NX-2-offset+1][jj  ][kk  ]    = E_Xdir_OLD[4+1][jj  ][kk  ]
                    + cnst * (    EB_WAVE[NX-4-offset+1][jj  ][kk  ] 
                              -E_Xdir_OLD[6+1                   ][jj  ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[NX-2-offset  ][jj+1][kk  ]    = E_Xdir_OLD[4  ][jj+1][kk  ]
                    + cnst * (    EB_WAVE[NX-4-offset  ][jj+1][kk  ] 
                              -E_Xdir_OLD[6                     ][jj+1][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[NX-2-offset  ][jj  ][kk+1]    = E_Xdir_OLD[4  ][jj  ][kk+1]
                    + cnst * (    EB_WAVE[NX-4-offset  ][jj  ][kk+1] 
                              -E_Xdir_OLD[6                     ][jj  ][kk+1] );
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
                EB_WAVE[ii+1][offset+0  ][kk  ] = E_Ydir_OLD[ii+1][2  ][kk  ]
                    + cnst * (    EB_WAVE[ii+1][offset+2  ][kk  ]
                              -E_Ydir_OLD[ii+1][0         ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[ii  ][offset+0+1][kk  ] = E_Ydir_OLD[ii  ][2+1][kk  ]
                    + cnst * (    EB_WAVE[ii  ][offset+2+1][kk  ]
                              -E_Ydir_OLD[ii  ][0+1       ][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[ii  ][offset+0  ][kk+1] = E_Ydir_OLD[ii  ][2  ][kk+1]
                    + cnst * (    EB_WAVE[ii  ][offset+2  ][kk+1]
                              -E_Ydir_OLD[ii  ][0         ][kk+1] );
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
                EB_WAVE[ii+1][NY-2-offset  ][kk  ] = E_Ydir_OLD[ii+1][4  ][kk  ]
                    + cnst * (    EB_WAVE[ii+1][NY-4-offset  ][kk  ]
                              -E_Ydir_OLD[ii+1][6                     ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[ii  ][NY-2-offset+1][kk  ] = E_Ydir_OLD[ii  ][4+1][kk  ]
                    + cnst * (    EB_WAVE[ii  ][NY-4-offset+1][kk  ]
                              -E_Ydir_OLD[ii  ][6+1                   ][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[ii  ][NY-2-offset  ][kk+1] = E_Ydir_OLD[ii  ][4  ][kk+1]
                    + cnst * (    EB_WAVE[ii  ][NY-4-offset  ][kk+1]
                              -E_Ydir_OLD[ii  ][6                     ][kk+1] );
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
                EB_WAVE[ii+1][jj  ][offset+0]   = E_Zdir_OLD[ii+1][jj  ][2  ]
                    + cnst * (    EB_WAVE[ii+1][jj  ][offset+2  ]
                              -E_Zdir_OLD[ii+1][jj  ][0  ]        );
                // Ey: even-odd-even
                EB_WAVE[ii  ][jj+1][offset+0]   = E_Zdir_OLD[ii  ][jj+1][2  ]
                    + cnst * (    EB_WAVE[ii  ][jj+1][offset+2  ]
                              -E_Zdir_OLD[ii  ][jj+1][0  ]        );
                // Ez: even-even-odd
                EB_WAVE[ii  ][jj  ][offset+0+1] = E_Zdir_OLD[ii  ][jj  ][2+1]
                    + cnst * (    EB_WAVE[ii  ][jj  ][offset+2+1]
                              -E_Zdir_OLD[ii  ][jj  ][0+1]        );
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
                EB_WAVE[ii+1][jj  ][NZ-2-offset  ]    = E_Zdir_OLD[ii+1][jj  ][4  ]
                    + cnst * (    EB_WAVE[ii+1][jj  ][NZ-4-offset  ]
                              -E_Zdir_OLD[ii+1][jj  ][6  ]     );
                // Ey: even-odd-even
                EB_WAVE[ii  ][jj+1][NZ-2-offset  ]    = E_Zdir_OLD[ii  ][jj+1][4  ]
                    + cnst * (    EB_WAVE[ii  ][jj+1][NZ-4-offset  ]
                              -E_Zdir_OLD[ii  ][jj+1][6  ]     );
                // Ez: even-even-odd
                EB_WAVE[ii  ][jj  ][NZ-2-offset+1]    = E_Zdir_OLD[ii  ][jj  ][4+1]
                    + cnst * (    EB_WAVE[ii  ][jj  ][NZ-4-offset+1]
                              -E_Zdir_OLD[ii  ][jj  ][6+1]     );
            }
        }
    }

    return EXIT_SUCCESS;

} //}}}


int abc_Mur_1st_ref( gridConfiguration *gridCfg,
                     double EB_WAVE[NX][NY][NZ_REF] ) {
//{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk,
        offset;             // refers to EB_WAVE only

    double
        cnst;

    cnst    = (DT-DX)/(DT+DX);
    offset  = 2;

    // absorber into x-direction
#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj<NY-2 ; jj+=2) {
        for (kk=2 ; kk<NZ_REF-2 ; kk+=2) {
            // absorber at x=0 grid boundary
            // Ex: odd-even-even
            EB_WAVE[offset+0+1][jj  ][kk  ] = E_Xdir_OLD_ref[2+1][jj  ][kk  ]
                + cnst * (    EB_WAVE[offset+2+1][jj  ][kk  ] 
                          -E_Xdir_OLD_ref[0+1       ][jj  ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[offset+0  ][jj+1][kk  ] = E_Xdir_OLD_ref[2  ][jj+1][kk  ]
                + cnst * (    EB_WAVE[offset+2  ][jj+1][kk  ] 
                          -E_Xdir_OLD_ref[0         ][jj+1][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[offset+0  ][jj  ][kk+1] = E_Xdir_OLD_ref[2  ][jj  ][kk+1]
                + cnst * (    EB_WAVE[offset+2  ][jj  ][kk+1] 
                          -E_Xdir_OLD_ref[0         ][jj  ][kk+1] );
            // absorber at x=NX grid boundary
            // Ex: odd-even-even
            EB_WAVE[NX-2-offset+1][jj  ][kk  ]    = E_Xdir_OLD_ref[4+1][jj  ][kk  ]
                + cnst * (    EB_WAVE[NX-4-offset+1][jj  ][kk  ] 
                          -E_Xdir_OLD_ref[6+1                   ][jj  ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[NX-2-offset  ][jj+1][kk  ]    = E_Xdir_OLD_ref[4  ][jj+1][kk  ]
                + cnst * (    EB_WAVE[NX-4-offset  ][jj+1][kk  ] 
                          -E_Xdir_OLD_ref[6                     ][jj+1][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[NX-2-offset  ][jj  ][kk+1]    = E_Xdir_OLD_ref[4  ][jj  ][kk+1]
                + cnst * (    EB_WAVE[NX-4-offset  ][jj  ][kk+1] 
                          -E_Xdir_OLD_ref[6                     ][jj  ][kk+1] );
        }
    }

    // absorber into y-direction
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (kk=2 ; kk<NZ_REF-2 ; kk+=2) {
            // absorber at y=0 grid boundary
            // Ex: odd-even-even
            EB_WAVE[ii+1][offset+0  ][kk  ] = E_Ydir_OLD_ref[ii+1][2  ][kk  ]
                + cnst * (    EB_WAVE[ii+1][offset+2  ][kk  ]
                          -E_Ydir_OLD_ref[ii+1][0         ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[ii  ][offset+0+1][kk  ] = E_Ydir_OLD_ref[ii  ][2+1][kk  ]
                + cnst * (    EB_WAVE[ii  ][offset+2+1][kk  ]
                          -E_Ydir_OLD_ref[ii  ][0+1       ][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[ii  ][offset+0  ][kk+1] = E_Ydir_OLD_ref[ii  ][2  ][kk+1]
                + cnst * (    EB_WAVE[ii  ][offset+2  ][kk+1]
                          -E_Ydir_OLD_ref[ii  ][0         ][kk+1] );
            // absorber at y=NY grid boundary
            // Ex: odd-even-even
            EB_WAVE[ii+1][NY-2-offset  ][kk  ] = E_Ydir_OLD_ref[ii+1][4  ][kk  ]
                + cnst * (    EB_WAVE[ii+1][NY-4-offset  ][kk  ]
                          -E_Ydir_OLD_ref[ii+1][6                     ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[ii  ][NY-2-offset+1][kk  ] = E_Ydir_OLD_ref[ii  ][4+1][kk  ]
                + cnst * (    EB_WAVE[ii  ][NY-4-offset+1][kk  ]
                          -E_Ydir_OLD_ref[ii  ][6+1                   ][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[ii  ][NY-2-offset  ][kk+1] = E_Ydir_OLD_ref[ii  ][4  ][kk+1]
                + cnst * (    EB_WAVE[ii  ][NY-4-offset  ][kk+1]
                          -E_Ydir_OLD_ref[ii  ][6                     ][kk+1] );
        }
    }

    // absorber into z-direction
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii<NX-2 ; ii+=2) {
        for (jj=2 ; jj<NY-2 ; jj+=2) {
            // absorber at z=0 grid boundary
            // Ex: odd-even-even
            EB_WAVE[ii+1][jj  ][offset+0]   = E_Zdir_OLD_ref[ii+1][jj  ][2  ]
                + cnst * (    EB_WAVE[ii+1][jj  ][offset+2  ]
                          -E_Zdir_OLD_ref[ii+1][jj  ][0  ]        );
            // Ey: even-odd-even
            EB_WAVE[ii  ][jj+1][offset+0]   = E_Zdir_OLD_ref[ii  ][jj+1][2  ]
                + cnst * (    EB_WAVE[ii  ][jj+1][offset+2  ]
                          -E_Zdir_OLD_ref[ii  ][jj+1][0  ]        );
            // Ez: even-even-odd
            EB_WAVE[ii  ][jj  ][offset+0+1] = E_Zdir_OLD_ref[ii  ][jj  ][2+1]
                + cnst * (    EB_WAVE[ii  ][jj  ][offset+2+1]
                          -E_Zdir_OLD_ref[ii  ][jj  ][0+1]        );
            // absorber at z=NZ grid boundary
            // Ex: odd-even-even
            EB_WAVE[ii+1][jj  ][NZ_REF-2-offset  ]    = E_Zdir_OLD_ref[ii+1][jj  ][4  ]
                + cnst * (    EB_WAVE[ii+1][jj  ][NZ_REF-4-offset  ]
                          -E_Zdir_OLD_ref[ii+1][jj  ][6  ]     );
            // Ey: even-odd-even
            EB_WAVE[ii  ][jj+1][NZ_REF-2-offset  ]    = E_Zdir_OLD_ref[ii  ][jj+1][4  ]
                + cnst * (    EB_WAVE[ii  ][jj+1][NZ_REF-4-offset  ]
                          -E_Zdir_OLD_ref[ii  ][jj+1][6  ]     );
            // Ez: even-even-odd
            EB_WAVE[ii  ][jj  ][NZ_REF-2-offset+1]    = E_Zdir_OLD_ref[ii  ][jj  ][4+1]
                + cnst * (    EB_WAVE[ii  ][jj  ][NZ_REF-4-offset+1]
                          -E_Zdir_OLD_ref[ii  ][jj  ][6+1]     );
        }
    }

    return EXIT_SUCCESS;

} //}}}

/*UPML functions*/
double sigma(int pml_size, double nn, int m, double ds){

    double sig, sig_max, R_0;
    
    if( pml_size == 20 ) R_0 = pow(10,-12);
    if( pml_size == 10 ) R_0 = pow(10,-6);
    if( pml_size == 5 ) R_0 = pow(10,-3);

    sig_max = -(m+1)*log( R_0 )/(2*ds*pml_size);
    sig = pow( (nn) /(pml_size), m) * sig_max;

    return sig;  
}

void init_UPML_parameters(   gridConfiguration *gridCfg, boundaryVariables *boundaryV){

    int ii, jj, kk, kkr, count;
    double sig, kx, ky, kz;
    
    kx = 1;
    ky = 1;
    kz = 1;

    count = (d_absorb-2)/2;
    for ( ii=2 ; ii < NX-2 ; ii+=2 ) {
        if(ii < d_absorb ){

            sig = sigma( (d_absorb-2)/2, count, 4, DX );
            F1x(ii/2) = (2*kx) - (sig*DT);
            F2x(ii/2) = (2*kx) + (sig*DT);
            Cx(ii/2) = F1x(ii/2)/F2x(ii/2);
            count -= 1;

        }else if( ii > NX - d_absorb - 2 ){

            count += 1;
            sig = sigma( (d_absorb-2)/2, count, 4, DX);
            F1x(ii/2) = (2*kx) - (sig*DT);
            F2x(ii/2) = (2*kx) + (sig*DT);
            Cx(ii/2) = F1x(ii/2)/F2x(ii/2);  

        }else{

            sig = 0;
            F1x(ii/2) = (2*kx) - (sig*DT);
            F2x(ii/2) = (2*kx) + (sig*DT);
            Cx(ii/2) = F1x(ii/2)/F2x(ii/2);
        }
        //printf("Cx(%d) = %.5f, count = %d \n", ii, Cx(ii/2), count );
    }
    
    count = (d_absorb-2)/2;
    for ( jj=2 ; jj < NY-2 ; jj+=2 ) {
        if(jj < d_absorb ){

            sig = sigma( (d_absorb-2)/2, count, 4, DX );
            F1y(jj/2) = (2*ky) - (sig*DT);
            F2y(jj/2) = (2*ky) + (sig*DT);
            Cy(jj/2) = F1y(jj/2)/F2y(jj/2);
            count -= 1;

        }else if( jj > NY - d_absorb - 2){

            count += 1;
            sig = sigma( (d_absorb-2)/2, count, 4, DX );
            F1y(jj/2) = (2*ky) - (sig*DT);
            F2y(jj/2) = (2*ky) + (sig*DT);
            Cy(jj/2) = F1y(jj/2)/F2y(jj/2); 

        }else{

            sig = 0;
            F1y(jj/2) = (2*ky) - (sig*DT);
            F2y(jj/2) = (2*ky) + (sig*DT);
            Cy(jj/2) = F1y(jj/2)/F2y(jj/2);
        }
        //printf("Cy(%d) = %.5f, count = %d \n", jj, Cy(jj/2), count );
    }
    
    count = (d_absorb-2)/2;
    for ( kk=2 ; kk < NZ-2 ; kk+=2 ) {
        if(kk < d_absorb){

            sig = sigma( (d_absorb-2)/2, count, 4, DX );
            F1z(kk/2) = (2*kz) - (sig*DT);
            F2z(kk/2) = (2*kz) + (sig*DT);
            Cz(kk/2) = F1z(kk/2)/F2z(kk/2);
            count -= 1;

        }else if( kk > NZ - d_absorb - 2){ 

            count += 1;
            sig = sigma( (d_absorb-2)/2 , count, 4, DX);
            F1z(kk/2) = (2*kz) - (sig*DT);
            F2z(kk/2) = (2*kz) + (sig*DT);
            Cz(kk/2) = F1z(kk/2)/F2z(kk/2); 

        }else{

            sig = 0;
            F1z(kk/2) = (2*kz) - (sig*DT);
            F2z(kk/2) = (2*kz) + (sig*DT);
            Cz(kk/2) = F1z(kk/2)/F2z(kk/2);
        }
        //printf("Cz(%d) = %.5f \n", kk, Cz(kk/2) );
    }

    count = (d_absorb-2)/2;
    for ( kkr=2 ; kkr < NZ_REF-2 ; kkr+=2 ) {
        if(kkr < d_absorb){

            sig = sigma( (d_absorb-2)/2, count, 4, DX );
            F1zr(kkr) = (2*kz) - (sig*DT);
            F2zr(kkr) = (2*kz) + (sig*DT);
            Czr(kkr) = F1zr(kkr)/F2zr(kkr);
            count -= 1;

        }else if( kkr > NZ_REF - d_absorb - 2 ){ 
                    
            count += 1;
            sig = sigma( (d_absorb-2)/2, count, 4, DX);
            F1zr(kkr) = (2*kz) - (sig*DT);
            F2zr(kkr) = (2*kz) + (sig*DT);
            Czr(kkr) = F1zr(kkr)/F2zr(kkr); 

        }else{
            sig = 0;
            F1zr(kkr) = (2*kz) - (sig*DT);
            F2zr(kkr) = (2*kz) + (sig*DT);
            Czr(kkr) = F1zr(kkr)/F2zr(kkr);
        }
        //printf("C(%d) = %.5f \n", kkr, Czr(kkr) );
    }

    //exit(-1);
}


int free_boundary(gridConfiguration *gridCfg){

    if(sel_boundary == 2){

        free(E_Xdir_OLD);

        printf("Mur boundary memory allocated free. \n");
    }
    

    return EXIT_SUCCESS;
}