#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "focal-struct.h"
#include "alloc-memory.h"
#include "macros-boundary.h"

int set_densityInAbsorber_v2( gridConfiguration *gridCfg,
                              char absorber[], 
                              double n_e[gridCfg->Nx/2][gridCfg->Ny/2][gridCfg->Nz/2] );

int apply_absorber( gridConfiguration *gridCfg, 
                    double eco, 
                    double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] );

int apply_absorber_ref( gridConfiguration *gridCfg, 
                        double eco, 
                        double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref] );

int apply_absorber_v2( size_t N_x, size_t N_y, size_t N_z, int d_absorb, double eco, 
                       char absorber[],
                       double EB_WAVE[N_x][N_y][N_z] );

int apply_numerical_viscosity( gridConfiguration *gridCfg,
                               double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] );

int abc_Mur_saveOldE_xdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[8][gridCfg->Ny][gridCfg->Nz] );
int abc_Mur_saveOldE_ydir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[gridCfg->Nx][8][gridCfg->Nz] );    // was [Ny,8,Nz] until 2022-01-22 (which is wrong)
int abc_Mur_saveOldE_zdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[gridCfg->Nx][gridCfg->Ny][8] );

int abc_Mur_saveOldEref_xdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[8][gridCfg->Ny][gridCfg->Nz] );
int abc_Mur_saveOldEref_ydir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[gridCfg->Nx][8][gridCfg->Nz] );
int abc_Mur_saveOldEref_zdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[gridCfg->Nx][gridCfg->Ny][8] );

int abc_Mur_1st( gridConfiguration *gridCfg, 
                 char absorber[],
                 double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                 double E_old_xdir[8][gridCfg->Ny][gridCfg->Nz], 
                 double E_old_ydir[gridCfg->Nx][8][gridCfg->Nz], 
                 double E_old_zdir[gridCfg->Nx][gridCfg->Ny][8] );
int abc_Mur_1st_ref( gridConfiguration *gridCfg, 
                     double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                     double E_old_xdir[8][gridCfg->Ny][gridCfg->Nz_ref], 
                     double E_old_ydir[gridCfg->Nx][8][gridCfg->Nz_ref], 
                     double E_old_zdir[gridCfg->Nx][gridCfg->Ny][8] );

//void initialize_split_PML(pmlBoundary *PML, gridConfiguration *gridCfg, int pml_size);

//void apply_split_PML(pmlBoundary *PML, gridConfiguration *gridCfg, int pml_size, double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz]);

//void split_PML_parameter(pmlBoundary *PML, gridConfiguration *gridCfg, int pml_size);

void setBoundary(gridConfiguration *gridCfg, namePath *pathFile, abcBoundary *ABC);

void computeBoundary(   gridConfiguration *gridCfg, 
                        namePath *pathFile,
                        abcBoundary *ABC,
                        double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], int t_int,
                        double EB_WAVE_ref [gridCfg->Nx][gridCfg->Ny][gridCfg->Nz]);

#endif