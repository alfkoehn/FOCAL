#ifndef ALLOC_MEMORY
#define ALLOC_MEMORY

#include <stdio.h>

#define ALLOC_1D(PNTR, NUM, TYPE)                                 \
    PNTR = (TYPE *)calloc( NUM, sizeof(TYPE));                    \
    if (!PNTR){                                                   \
        perror("ALLOC_1D");                                       \
        fprintf(stderr,                                           \
          "Allocation failed for " #PNTR ". Terminating...\n");   \
        exit(-1);                                                 \
    }

#define ALLOC_2D(PNTR, NUMX, NUMY, TYPE)                          \
    PNTR = (TYPE *)calloc( (NUMX)*(NUMY), sizeof(TYPE));          \
    if (!PNTR){                                                   \
        perror("ALLOC_2D");                                       \
        fprintf(stderr,                                           \
          "Allocation failed for " #PNTR ". Terminating...\n");   \
        exit(-1);                                                 \
    }

#define ALLOC_3D(PNTR, NUMX, NUMY, NUMZ, TYPE)                    \
    PNTR = (TYPE *)calloc( (NUMX)*(NUMY)*(NUMZ), sizeof(TYPE));   \
    if (!PNTR){                                                   \
        perror("ALLOC_3D");                                       \
        fprintf(stderr,                                           \
          "Allocation failed for " #PNTR ". Terminating...\n");   \
        exit(-1);                                                 \
    }


#endif