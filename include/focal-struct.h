#ifndef FOCAL_STRUCT_H
#define FOCAL_STRUCT_H

// define structures
typedef struct gridConfiguration {
    int
        Nx, Ny, Nz,     // maybe size_t would be better
        Nz_ref,
        d_absorb,
        t_end,
        ne_profile, B0_profile;
    double
        period,
        dx,dt,
        ne_max, B0_value;
} gridConfiguration;
typedef struct beamConfiguration {
    int
        exc_signal,
        ant_x, ant_y, ant_z,
        rampUpMethod;
    double
        antAngle_zy, antAngle_zx,
        ant_w0x, ant_w0y,
        z2waist;
} beamConfiguration;
typedef struct namePath{
    const char 	*foldername, *projectPath,
    		    *file_hdf5, *file_trace;
    int boundary;
} namePath;
typedef struct antennaDetector{
    int antDetect_1D, 
        detAnt_01_zpos, detAnt_02_zpos,
        detAnt_03_zpos,detAnt_04_zpos,
        detAnt_01_ypos;
} antennaDetector;
typedef struct abcBoundary{
    double eco;
} abcBoundary;


#endif