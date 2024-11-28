// header guard first to prevent multiple declarations
#ifndef FOCAL_H
#define FOCAL_H

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "focal-struct.h"
#include "macros-grid.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

void advance_fields(    gridConfiguration *gridCfg, 
                        double EB_WAVE[NX][NY][NZ], 
                        double EB_WAVE_ref[NX][NY][NZ_REF],
                        double J_B0[NX][NY][NZ],
                        double n_e[NX/2][NY/2][NZ/2] );


int advance_J( gridConfiguration *gridCfg, 
               double EB_WAVE[NX][NY][NZ], 
               double J_B0[NX][NY][NZ],
               double n_e[NX/2][NY/2][NZ/2] ); 

int advance_B( gridConfiguration *gridCfg, 
               double EB_WAVE[NX][NY][NZ] );

int advance_B_ref( gridConfiguration *gridCfg, 
                   double EB_WAVE[NX][NY][NZ_REF] );

int advance_E( gridConfiguration *gridCfg, 
               double EB_WAVE[NX][NY][NZ], 
               double J_B0[NX][NY][NZ] );

int advance_E_ref( gridConfiguration *gridCfg, 
                   double EB_WAVE[NX][NY][NZ_REF] ); 

int set_densityInAbsorber_v2( gridConfiguration *gridCfg,
                              char absorber[], 
                              double n_e[NX/2][NY/2][NZ/2] );

int apply_numerical_viscosity( gridConfiguration *gridCfg,
                               double EB_WAVE[NX][NY][NZ] );

int set2zero_1D( size_t N_x, double arr_1D[N_x] );
int set2zero_3D( size_t N_x, size_t N_y, size_t N_z, double arr_3D[N_x][N_y][N_z] );

int advance_J_PML( gridConfiguration *gridCfg, 
               double EB_WAVE[NX][NY][NZ], 
               double J_B0[NX][NY][NZ],
               double n_e[NX/2][NY/2][NZ/2] );

int advance_B_PML(  gridConfiguration *gridCfg, 
                    double EB_WAVE[NX][NY][NZ] );

int advance_Bref_PML(   gridConfiguration *gridCfg, 
                        double EB_WAVE[NX][NY][NZ_REF] );

int advance_E_PML(  gridConfiguration *gridCfg, 
                    double EB_WAVE[NX][NY][NZ], 
                    double J_B0[NX][NY][NZ] );

int advance_Eref_PML(   gridConfiguration *gridCfg, 
                        double EB_WAVE[NX][NY][NZ_REF] );

#endif  // FOCAL_H

