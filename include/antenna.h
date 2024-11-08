#ifndef ANTENNA_H
#define ANTANNA_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "focal-struct.h"
#include "macros-grid.h"
//#include "../include/focal.h"
#include "focal.h"

int make_antenna_profile( gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg,
                          double antField_xy[NX/2][Ny/2], double antPhaseTerms[NX/2][Ny/2] );

int add_source( gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg, 
                int t_int, double omega_t, 
                double antField_xy[NX/2][Ny/2], 
                double antPhaseTerms[NX/2][Ny/2],
                double EB_WAVE[NX][Ny][Nz] );

int add_source_ref( gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg, 
                    int t_int, double omega_t, 
                    double antField_xy[NX/2][Ny/2], 
                    double antPhaseTerms[NX/2][Ny/2],
                    double EB_WAVE[NX][Ny][Nz_ref] );

double antenna_field_rampup( int RampUpMethod, double Period, int t_int );

double antenna_calcHansenExEy_O( double theta_rad, double Y );

#endif  // ANTENNA_H
