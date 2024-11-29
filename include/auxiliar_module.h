#ifndef AUXILIAR_MODULE_H
#define AUXILIAR_MODULE_H

#include <stdlib.h>
#include <stdio.h>

// Function to allocate memory for a 3D array with given dimensions
double ***allocate3DArray(int N_x, int N_y, int N_z);
void free3DArray(double ***array, int N_x, int N_y);

// Function to allocate memory for a 2D array with given dimensions
double **allocate2DArray(int N_x, int N_y);
void free2DArray(double **array, int N_x);

#endif