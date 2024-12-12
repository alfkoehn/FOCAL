#ifndef POWER_CALC_H
#define POWER_CALC_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "focal-struct.h"
#include "macros-grid.h"
#include "grid_io.h"
#include "auxiliar_module.h"

//Initialization of the power related values
int init_powerValues(   gridConfiguration *gridCfg,
                        powerValues *powerVal );

//function called in MAIN to compute power
void compute_power( gridConfiguration *gridCfg,
                    beamAntennaConfiguration *beamCfg,
                    powerValues *powerVal,
                    int t_int, 
                    double EB_WAVE[NX][NY][NZ], 
                    double EB_WAVE_ref[NX][NY][NZ_REF] );

//calls power calculation
int calculate_power(gridConfiguration *gridCfg, 
                    powerValues *powerVal,
                    int t_int, 
                    double EB_WAVE[NX][NY][NZ], 
                    double EB_WAVE_ref[NX][NY][NZ_REF] );

//Store power values in Timetraces
int power_toTimetraces( gridConfiguration *gridCfg, 
                        beamAntennaConfiguration *beamCfg, 
                        powerValues *powerVal,
                        int t_int );

//Functions to compute power
double calc_poynt_4( gridConfiguration *gridCfg, 
                     powerValues *powerVal, char absorber[],
                     double EB_WAVE[NX][NY][NZ], 
                     double EB_WAVE_ref[NX][NY][NZ_REF] );
double calc_poynt_5( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     powerValues *powerVal, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] );
double calc_poynt_6( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     powerValues *powerVal, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] );
double calc_poynt_7( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     powerValues *powerVal, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] );

//Print time traces to console and .file
int write_timetraces(   gridConfiguration *gridCfg,
                        saveData *saveDCfg );
//Print full timetraces pointer to console
int writeConsole_timetraces( int T_end, double Period );

#endif  // POWER_CALC_H
