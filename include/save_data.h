#ifndef SAVE_DATA_H
#define SAVE_DATA_H

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#include "focal-struct.h"
#include "macros-grid.h"

void create_folder(saveData *saveDCfg);
void simulation_folder(const char *path);
void data_folder(const char *path, const char *folder_name);
void copyJSON(const char *path, const char *folder_name);

#endif