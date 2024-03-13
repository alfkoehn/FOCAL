#ifndef GRID_IO_H
#define GRID_IO_H

int writeTimetraces2ascii( int dim0, int dim1, int t_end, double period, 
                           char filename[], double timetraces[dim0][dim1] );

//#ifdef HDF5
int writeMyHDF_v4( int dim0, int dim1, int dim2, char filename[], char dataset[], double array_3D[dim0][dim1][dim2] );
int writeConfig2HDF( gridConfiguration *gridCfg, beamConfiguration *beamCfg, char filename[] );
int readMyHDF( int dim0, int dim1, int dim2, char filename[], char dataset[], double array_3D[dim0][dim1][dim2]);
//#endif

//#ifdef DETECTOR_ANTENNA_1D
int detAnt1D_storeValues( gridConfiguration *gridCfg,
                          size_t detAnt_ypos, size_t detAnt_zpos,
                          int tt, 
                          double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                          double detAnt_fields[gridCfg->Nx/2][5] );
//#endif

//#if defined(HDF5) && defined(DETECTOR_ANTENNA_1D)
int detAnt1D_write2hdf5( int N_x, 
                         char filename[], char detAnt_groupName[], 
                         size_t detAnt_ypos, size_t detAnt_zpos,
                         double detAnt_fields[N_x/2][5] );
//#endif



#endif  // GRID_IO_H
