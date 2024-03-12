#ifndef GRID_IO_H
#define GRID_IO_H

int writeTimetraces2ascii( int dim0, int dim1, int t_end, double period, 
                           char filename[], double timetraces[dim0][dim1] );

#ifdef HDF5
int writeMyHDF_v4( int dim0, int dim1, int dim2, char filename[], char dataset[], double array_3D[dim0][dim1][dim2] );
int writeConfig2HDF( gridConfiguration *gridCfg, beamConfiguration *beamCfg, char filename[] );
int readMyHDF( int dim0, int dim1, int dim2, char filename[], char dataset[], double array_3D[dim0][dim1][dim2]);
#endif

#endif  // GRID_IO_H
