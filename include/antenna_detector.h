#ifndef ANTENNA_DETECTOR_H
#define ANTENNA_DETECTOR_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "focal-struct.h"
#include "macros-grid.h"
#include "grid_io.h"
#include "auxiliar_module.h"

int init_antennaDetect( gridConfiguration *gridCfg,
                        beamAntennaConfiguration *beamCfg,
                        antennaDetector *antDetect );

int free_antDetect( gridConfiguration *gridCfg,
                    antennaDetector *antDetect );

int print_antennaDetec( antennaDetector *antDetect );

int control_antennaDetect(  gridConfiguration *gridCfg,
                            antennaDetector *antDetect,
                            int t_int,
                            double EB_WAVE[NX][NY][NZ] );

int detAnt1D_storeValues(   gridConfiguration *gridCfg, 
                            size_t detAnt_ypos, size_t detAnt_zpos,
                            int tt, 
                            double EB_WAVE[NX][NY][NZ], 
                            double **detAnt_fields );

void save_AntDetect(    gridConfiguration *gridCfg, saveData *saveDCfg,
                        antennaDetector *antDetect );

#endif