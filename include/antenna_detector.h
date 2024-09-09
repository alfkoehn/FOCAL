#ifndef ANTENNA_DETECTOR_H
#define ANTENNA_DETECTOR_H

#include <math.h>
#include <stdio.h>

#include "focal-struct.h"

void initialize_antDetect(antennaDetector *ant_Detect, gridConfiguration *gridCfg, beamConfiguration *beamCfg);
void print_antDetect_info(antennaDetector *ant_Detect);

#endif