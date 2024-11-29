#include "auxiliar_module.h"

// Function to allocate memory for a 3D array with given dimensions
double ***allocate3DArray(int N_x, int N_y, int N_z) {

    double ***array = (double ***)calloc( N_x, sizeof(double **));

#pragma omp parallel for
    for (int i = 0; i < N_x; i++) {
        array[i] = (double **)calloc( N_y, sizeof(double *));
        for (int j = 0; j < N_y; j++) {
            array[i][j] = (double *)calloc( N_z, sizeof(double));
        }
    }
    
    return array;
}

void free3DArray(double ***array, int N_x, int N_y){

    int i, j;
#pragma omp parallel for private(j)
    for( i = 0 ; i<N_x ; i+=1 ){
        for( j = 0 ; j<N_y ; j+=1 ){
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}

// Function to allocate memory for a 2D array with given dimensions
double **allocate2DArray(int N_x, int N_y) {

    double **array = (double **)calloc( N_x, sizeof(double *));

#pragma omp parallel for
    for (int i = 0; i < N_x; i++) {
        array[i] = (double *)calloc( N_y, sizeof(double ));
    }
    
    return array;
}

void free2DArray(double **array, int N_x){

#pragma omp parallel for
    for(int i = 0 ; i<N_x ; i+=1 ){
        free(array[i]);
    }
    free(array);
}

int set2zero3DArray( double ***arr_3D, size_t N_x, size_t N_y, size_t N_z  ){
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
