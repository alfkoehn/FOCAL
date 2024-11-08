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

#define ABSORBER_DAMPING(eco,damp) (1.-eco*damp*damp)


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

int apply_absorber( gridConfiguration *gridCfg, 
                    double eco, 
                    double EB_WAVE[NX][NY][NZ] );

int apply_absorber_ref( gridConfiguration *gridCfg, 
                        double eco, 
                        double EB_WAVE[NX][NY][NZ_REF] );

int apply_absorber_v2( size_t N_x, size_t N_y, size_t N_z, int D_absorb, double eco, 
                       char absorber[],
                       double EB_WAVE[N_x][N_y][N_z] );

int apply_numerical_viscosity( gridConfiguration *gridCfg,
                               double EB_WAVE[NX][NY][NZ] );

int abc_Mur_saveOldE_xdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[NX][NY][NZ], 
                           double E_old[8][NY][NZ] );
int abc_Mur_saveOldE_ydir( gridConfiguration *gridCfg, 
                           double EB_WAVE[NX][NY][NZ], 
                           double E_old[NX][8][NZ] );    // was [NY,8,NZ] until 2022-01-22 (which is wrong)
int abc_Mur_saveOldE_zdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[NX][NY][NZ], 
                           double E_old[NX][NY][8] );

int abc_Mur_saveOldEref_xdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[NX][NY][NZ_REF], 
                              double E_old[8][NY][NZ] );
int abc_Mur_saveOldEref_ydir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[NX][NY][NZ_REF], 
                              double E_old[NX][8][NZ] );
int abc_Mur_saveOldEref_zdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[NX][NY][NZ_REF], 
                              double E_old[NX][NY][8] );

int abc_Mur_1st( gridConfiguration *gridCfg, 
                 char absorber[],
                 double EB_WAVE[NX][NY][NZ], 
                 double E_old_xdir[8][NY][NZ], 
                 double E_old_ydir[NX][8][NZ], 
                 double E_old_zdir[NX][NY][8] );

int abc_Mur_1st_ref( gridConfiguration *gridCfg, 
                     double EB_WAVE[NX][NY][NZ_REF], 
                     double E_old_xdir[8][NY][NZ_REF], 
                     double E_old_ydir[NX][8][NZ_REF], 
                     double E_old_zdir[NX][NY][8] );

int set2zero_1D( size_t N_x, double arr_1D[N_x] );
int set2zero_3D( size_t N_x, size_t N_y, size_t N_z, double arr_3D[N_x][N_y][N_z] );

#endif  // FOCAL_H

