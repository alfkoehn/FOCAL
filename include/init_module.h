#ifndef INIT_MODULE_H
#define INIT_MODULE_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>

#include "cJSON.h"
#include "focal-struct.h"
#include "macros-grid.h"

void control_init(  gridConfiguration *gridCfg, 
                    beamAntennaConfiguration *beamCfg,
                    saveData *saveDCfg );

void grid_init( gridConfiguration *gridCfg, 
                beamAntennaConfiguration *beamCfg,
                saveData *saveDCfg );

/*Functions in charge of JSON reading*/
void write_JSON_toGrid( gridConfiguration *gridCfg, 
                        beamAntennaConfiguration *beamCfg,
                        saveData *saveDCfg );

char *read_json();

void print_systemConfiguration(gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg );

#endif