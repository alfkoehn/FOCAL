#ifndef BOUNDARY_MODULE_H
#define BOUNDARY_MODULE_H

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "focal-struct.h"
#include "macros-grid.h"
#include "alloc-memory.h"
#include "UPML_module.h"

#define ABSORBER_DAMPING(eco,damp) (1.-eco*damp*damp)

void init_boundary(gridConfiguration *gridCfg, boundaryVariables *boundaryV);

double ***allocateBoundaryArray(int N_x, int N_y, int N_z);

void advance_boundary(  gridConfiguration *gridCfg, boundaryVariables *boundaryV, 
                        double EB_WAVE[NX][NY][NZ], double EB_WAVE_ref[NX][NY][NZ_REF]);

/*ABC functions*/
int apply_absorber( gridConfiguration *gridCfg, 
                    boundaryVariables *boundaryV, 
                    double EB_WAVE[NX][NY][NZ] );

int apply_absorber_ref( gridConfiguration *gridCfg, 
                        boundaryVariables *boundaryV, 
                        double EB_WAVE[NX][NY][NZ_REF] );

int apply_absorber_v2(  size_t N_x, size_t N_y, size_t N_z, int D_absorb, 
                        boundaryVariables *boundaryV, 
                        char absorber[],
                        double EB_WAVE[N_x][N_y][N_z] );

/*Mur boundary functions*/
int abc_Mur_saveOldE_xdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[NX][NY][NZ] );

int abc_Mur_saveOldE_ydir( gridConfiguration *gridCfg, 
                           double EB_WAVE[NX][NY][NZ] );

int abc_Mur_saveOldE_zdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[NX][NY][NZ] );

int abc_Mur_saveOldEref_xdir(   gridConfiguration *gridCfg, 
                                double EB_WAVE_ref[NX][NY][NZ_REF] );

int abc_Mur_saveOldEref_ydir(   gridConfiguration *gridCfg, 
                                double EB_WAVE_ref[NX][NY][NZ_REF] );

int abc_Mur_saveOldEref_zdir(   gridConfiguration *gridCfg, 
                                double EB_WAVE_ref[NX][NY][NZ_REF] );

int abc_Mur_1st( gridConfiguration *gridCfg, 
                 char absorber[],
                 double EB_WAVE[NX][NY][NZ] );

int abc_Mur_1st_ref( gridConfiguration *gridCfg,
                     double EB_WAVE[NX][NY][NZ_REF] );

/*UPML functions*/
double sigma(int pml_size, double nn, int m, double ds);

void init_UPML_parameters(   gridConfiguration *gridCfg, boundaryVariables *boundaryV);

#endif