#ifndef SAVE_DATA_H
#define SAVE_DATA_H

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include "focal.h"

void simulation_folder(const char *path);
void data_folder(const char *path, const char *foldername);
void copyJSON(const char *path, const char *foldername);
void create_folder_path(namePath *pathFile);

#endif
