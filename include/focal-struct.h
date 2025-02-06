#ifndef FOCAL_STRUCT_H
#define FOCAL_STRUCT_H

// define structures
typedef struct gridConfiguration {      /*Saves the main physical parameters of the system*/
    int
        Nx, Ny, Nz,   
        Nz_ref,
        d_absorb,
        t_end,
        ne_profile, B0_profile,
        sel_boundary;
    double
        period,
        dx,dt,
        ne_0, B0_value;
} gridConfiguration;

typedef struct saveData{                /*Variables related to simulation data saving*/
    
    const char 	*foldername, *projectPath,
    		    *file_hdf5, *file_trace,
                *file_input;

    int     t_save, col_for_timetraces;
} saveData;

typedef struct boundaryVariables{            /*Store grid value for the boundary variables*/

    /*ABC boundary*/
    double eco,

    /*UPML boundary*/
    *F1x, *F1y, *F1z,
    *F2x, *F2y, *F2z,
    *Cx, *Cy, *Cz,
    *F1zr, *F2zr,
    *Czr;
} boundaryVariables;

typedef struct beamAntennaConfiguration {   /*Antenna configuration variables*/
    int
        T_wave,
        exc_signal,
        ant_x, ant_y, ant_z,
        rampUpMethod;
    double
        omega_t,
        antAngle_zy, antAngle_zx,
        ant_w0x, ant_w0y,
        z2waist,
        R_pola_xy,                          // ratio of axes in polarization ellipse, see Hansen1985PPCF, eq. (A6)
        Y_at_X1, k0Ln_at_X1, theta_at_X1;
} beamAntennaConfiguration;

typedef struct antennaDetector{
    int 
        antDetect_1D, 
        detAnt_01_zpos, detAnt_02_zpos,
        detAnt_03_zpos, detAnt_04_zpos,
        detAnt_01_ypos;
} antennaDetector;

typedef struct powerValues{
    int 
        pwr_dect;

    double
        power_abs_x1, power_abs_x2,
        power_abs_y1, power_abs_y2,
        power_abs_z1, power_abs_z2,
        power_abs_ref,

        poynt_x1, poynt_x2,
        poynt_y1, poynt_y2,
        poynt_z1, poynt_z2,
        poynt_z1_ref;    
} powerValues;

#endif
