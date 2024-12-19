#ifndef SAVE_DATA_H
#define SAVE_DATA_H

#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#include "focal-struct.h"
#include "macros-grid.h"
#include "focal.h"
#include "grid_io.h"
#include "auxiliar_module.h"

void create_folder(gridConfiguration *gridCfg, saveData *saveDCfg);
void simulation_folder(const char *path);
void data_folder(const char *path, const char *folder_name);
void copyJSON(const char *path, const char *folder_name);

/*Functions to save data in folders*/
void save_SimData(  gridConfiguration *gridCfg,
                    beamAntennaConfiguration *beamCfg, 
                    saveData *saveDCfg,
                    double n_e[NX/2][NY/2][NZ/2],
                    double J_B0[NX][NY][NZ] );

int save_data_toHDF5(   gridConfiguration *gridCfg,
                        beamAntennaConfiguration *beamCfg,
                        char filename_hdf5[],
                        double n_e[NX/2][NY/2][NZ/2],
                        double J_B0[NX][NY][NZ] );

int save_field_toHDF5(  gridConfiguration *gridCfg, 
                        saveData *saveDCfg, int t_int,
                        double EB_WAVE[NX][NY][NZ] );

#endif