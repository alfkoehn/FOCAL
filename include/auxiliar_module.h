#ifndef AUXILIAR_MODULE_H
#define AUXILIAR_MODULE_H

#include <stdlib.h>

// Function to allocate memory for a 3D array with given dimensions
double ***allocate3DArray(int N_x, int N_y, int N_z);
void free3DArray(double ***array, int N_x, int N_y);

#endif