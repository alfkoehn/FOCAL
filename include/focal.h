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
               double EB_WAVE[Nx][Ny][Nz], 
               double J_B0[Nx][Ny][Nz],
               double n_e[Nx/2][Ny/2][Nz/2] ); 

int advance_B( gridConfiguration *gridCfg, 
               double EB_WAVE[Nx][Ny][Nz] );

int advance_B_ref( gridConfiguration *gridCfg, 
                   double EB_WAVE[Nx][Ny][Nz_ref] );

int advance_E( gridConfiguration *gridCfg, 
               double EB_WAVE[Nx][Ny][Nz], 
               double J_B0[Nx][Ny][Nz] );

int advance_E_ref( gridConfiguration *gridCfg, 
                   double EB_WAVE[Nx][Ny][Nz_ref] ); 

int set_densityInAbsorber_v2( gridConfiguration *gridCfg,
                              char absorber[], 
                              double n_e[Nx/2][Ny/2][Nz/2] );

int apply_absorber( gridConfiguration *gridCfg, 
                    double eco, 
                    double EB_WAVE[Nx][Ny][Nz] );

int apply_absorber_ref( gridConfiguration *gridCfg, 
                        double eco, 
                        double EB_WAVE[Nx][Ny][Nz_ref] );

int apply_absorber_v2( size_t N_x, size_t N_y, size_t N_z, int D_absorb, double eco, 
                       char absorber[],
                       double EB_WAVE[N_x][N_y][N_z] );

int apply_numerical_viscosity( gridConfiguration *gridCfg,
                               double EB_WAVE[Nx][Ny][Nz] );

int abc_Mur_saveOldE_xdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[Nx][Ny][Nz], 
                           double E_old[8][Ny][Nz] );
int abc_Mur_saveOldE_ydir( gridConfiguration *gridCfg, 
                           double EB_WAVE[Nx][Ny][Nz], 
                           double E_old[Nx][8][Nz] );    // was [Ny,8,Nz] until 2022-01-22 (which is wrong)
int abc_Mur_saveOldE_zdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[Nx][Ny][Nz], 
                           double E_old[Nx][Ny][8] );

int abc_Mur_saveOldEref_xdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[Nx][Ny][Nz_ref], 
                              double E_old[8][Ny][Nz] );
int abc_Mur_saveOldEref_ydir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[Nx][Ny][Nz_ref], 
                              double E_old[Nx][8][Nz] );
int abc_Mur_saveOldEref_zdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[Nx][Ny][Nz_ref], 
                              double E_old[Nx][Ny][8] );

int abc_Mur_1st( gridConfiguration *gridCfg, 
                 char absorber[],
                 double EB_WAVE[Nx][Ny][Nz], 
                 double E_old_xdir[8][Ny][Nz], 
                 double E_old_ydir[Nx][8][Nz], 
                 double E_old_zdir[Nx][Ny][8] );

int abc_Mur_1st_ref( gridConfiguration *gridCfg, 
                     double EB_WAVE[Nx][Ny][Nz_ref], 
                     double E_old_xdir[8][Ny][Nz_ref], 
                     double E_old_ydir[Nx][8][Nz_ref], 
                     double E_old_zdir[Nx][Ny][8] );

int set2zero_1D( size_t N_x, double arr_1D[N_x] );
int set2zero_3D( size_t N_x, size_t N_y, size_t N_z, double arr_3D[N_x][N_y][N_z] );

#endif  // FOCAL_H

