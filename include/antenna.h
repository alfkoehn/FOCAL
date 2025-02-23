#ifndef ANTENNA_H
#define ANTANNA_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "focal-struct.h"
#include "macros-grid.h"
#include "auxiliar_module.h"

void init_antennaInjection( gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg );
int make_antenna_profile(   gridConfiguration *gridCfg, 
                            beamAntennaConfiguration *beamCfg );

void control_antennaInjection(  gridConfiguration *gridCfg, 
                                beamAntennaConfiguration *beamCfg,
                                int t_int,
                                double EB_WAVE[NX][NY][NZ],
                                double EB_WAVE_ref[NX][NY][NZ_REF] );

int add_source( gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg, 
                int t_int, 
                double EB_WAVE[NX][NY][NZ] );

int add_source_ref( gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg, 
                    int t_int, 
                    double EB_WAVE[NX][NY][NZ_REF] );

int add_source_sourceRef( gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg, 
                          int t_int, 
                          double EB_WAVE[NX][NY][NZ], double EB_WAVE_ref[NX][NY][NZ_REF] );

double antenna_field_rampup( int RampUpMethod, double period, int t_int );

double antenna_calcHansenExEy_O( double theta_rad, double Y );

#endif  // ANTENNA_H
