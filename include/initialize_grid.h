#ifndef INITIALIZE_GRID_H
#define INITIALIZE_GRID_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cJSON.h"
#include "focal.h"

//functions in initialize grid
void gridConfInit(gridConfiguration *gridCfg, namePath *pathFile, beamConfiguration *beamCfg, antennaDetector *ant_Detect);
char *read_json();
void write_JSON_onGrid(gridConfiguration *gridCfg, namePath *pathFile, beamConfiguration *beamCfg, antennaDetector *ant_Detect);

#endif