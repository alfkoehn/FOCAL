#ifndef POWER_CALC_H
#define POWER_CALC_H

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "focal-struct.h"
#include "focal.h"

double calc_poynt_4( gridConfiguration *gridCfg, 
                     int pwr_dect, char absorber[],
                     double EB_WAVE[NX][NY][NZ], 
                     double EB_WAVE_ref[NX][NY][NZ_REF] );
double calc_poynt_5( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] );
double calc_poynt_6( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] );
double calc_poynt_7( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] );

#endif  // POWER_CALC_H
