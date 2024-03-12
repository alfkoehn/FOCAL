/**
 * Author:      Alf Köhn-Seemann
 * Email:       koehn@igvp.uni-stuttgart.de
 * Copyright:   University of Stuttgart
 * 
 * This is a 3D FDTD code for simulating electromagnetic waves in cold 
 * magnetized plasmas.
 *
 * NOTE: This is an early version, including some obsolete function, those
 *       will be removed in near future.
 *       Furthermore, everything will be properly split into separate 
 *       libraries, allowing the usage of a nice make file.
 *
 * Initial release on github: 2022-03-31
 *
 **/

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <stdarg.h>
#include <getopt.h>
#include <sys/stat.h>
#include <stdbool.h>


// check if compiler understands OMP, if not, this file does probably not exist
#ifdef _OPENMP
    #include <omp.h>  
#endif

#define HDF5
#ifdef HDF5
    #include "hdf5.h"
#endif

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

#define ABSORBER_DAMPING(eco,damp) (1.-eco*damp*damp)

// setting boundary conditions, possible choices are
// 1: simple_abc
// 2: Mur
#define BOUNDARY 1

#define DETECTOR_ANTENNA_1D


#include "focal.h"
#include "antenna.h"
#include "grid_io.h"


// prototyping
int make_density_profile( gridConfiguration *gridCfg, 
                          double cntrl_para, 
                          double n_e[gridCfg->Nx/2][gridCfg->Ny/2][gridCfg->Nz/2] );
int set_densityInAbsorber_v2( gridConfiguration *gridCfg,
                              char absorber[], 
                              double n_e[gridCfg->Nx/2][gridCfg->Ny/2][gridCfg->Nz/2] );
int make_B0_profile( gridConfiguration *gridCfg,
                     double cntrl_para, 
                     double J_B0[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] );
int apply_absorber( gridConfiguration *gridCfg, 
                    double eco, 
                    double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] );
int apply_absorber_ref( gridConfiguration *gridCfg, 
                        double eco, 
                        double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref] );
int apply_absorber_v2( size_t N_x, size_t N_y, size_t N_z, int d_absorb, double eco, 
                       char absorber[],
                       double EB_WAVE[N_x][N_y][N_z] );
int apply_numerical_viscosity( gridConfiguration *gridCfg,
                               double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] );
int abc_Mur_saveOldE_xdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[8][gridCfg->Ny][gridCfg->Nz] );
int abc_Mur_saveOldE_ydir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[gridCfg->Nx][8][gridCfg->Nz] );    // was [Ny,8,Nz] until 2022-01-22 (which is wrong)
int abc_Mur_saveOldE_zdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[gridCfg->Nx][gridCfg->Ny][8] );
int abc_Mur_saveOldEref_xdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[8][gridCfg->Ny][gridCfg->Nz] );
int abc_Mur_saveOldEref_ydir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[gridCfg->Nx][8][gridCfg->Nz] );
int abc_Mur_saveOldEref_zdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[gridCfg->Nx][gridCfg->Ny][8] );
int abc_Mur_1st( gridConfiguration *gridCfg, 
                 char absorber[],
                 double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                 double E_old_xdir[8][gridCfg->Ny][gridCfg->Nz], 
                 double E_old_ydir[gridCfg->Nx][8][gridCfg->Nz], 
                 double E_old_zdir[gridCfg->Nx][gridCfg->Ny][8] );
int abc_Mur_1st_ref( gridConfiguration *gridCfg, 
                     double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                     double E_old_xdir[8][gridCfg->Ny][gridCfg->Nz_ref], 
                     double E_old_ydir[gridCfg->Nx][8][gridCfg->Nz_ref], 
                     double E_old_zdir[gridCfg->Nx][gridCfg->Ny][8] );
double calc_poynt_4( gridConfiguration *gridCfg, 
                     int pwr_dect, char absorber[],
                     double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                     double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref] );
double calc_poynt_5( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] );
double calc_poynt_6( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] );
double calc_poynt_7( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] );
double calc_power_EE_1( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                       int pwr_dect, char absorber[],
                       double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] );
int set2zero_1D( size_t N_x, double arr_1D[N_x] );
int set2zero_3D( size_t N_x, size_t N_y, size_t N_z, double arr_3D[N_x][N_y][N_z] );
#ifdef DETECTOR_ANTENNA_1D
int detAnt1D_storeValues( gridConfiguration *gridCfg,
                          size_t detAnt_ypos, size_t detAnt_zpos,
                          int tt, 
                          double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                          double detAnt_fields[gridCfg->Nx/2][5] );
#endif
#if defined(HDF5) && defined(DETECTOR_ANTENNA_1D)
int detAnt1D_write2hdf5( int N_x, 
                         char filename[], char detAnt_groupName[], 
                         size_t detAnt_ypos, size_t detAnt_zpos,
                         double detAnt_fields[N_x/2][5] );
#endif


int main( int argc, char *argv[] ) {
//{{{

    struct gridConfiguration gridCfg;
    struct beamConfiguration beamCfg;

    int
        ii,jj,kk,
        t_int, T_wave, 

        scale,

#ifdef _OPENMP
        n_threads,                          // number of threads that will be used (OpenMP)
#endif
        pwr_dect,

#ifdef DETECTOR_ANTENNA_1D
        detAnt_01_zpos,
        detAnt_02_zpos,
        detAnt_03_zpos,
        detAnt_04_zpos,
        detAnt_01_ypos,
#endif

        len_str,                            // return value of snprintf
        opt_ret;                            // return value of getopt (reading input parameter)

    double
#if BOUNDARY == 1
        eco,
#endif

        //ant_phase, 

        poynt_x1, poynt_x2,
        poynt_y1, poynt_y2,
        poynt_z1, poynt_z2,
        poynt_z1_ref,

        power_abs_x1, power_abs_x2,
        power_abs_y1, power_abs_y2,
        power_abs_z1, power_abs_z2,
        power_abs_ref,

        /*
        power_EE_x1, power_EE_x2, 
        power_EE_y1, power_EE_y2, 
        power_EE_z1, power_EE_z2, 
        power_EE_ref, 
        */

        omega_t;

    char
        dSet_name[PATH_MAX],
        filename_hdf5[PATH_MAX];            // filename of hdf5 file for output

    bool
        angle_zx_set,                       // is antAngle_zx set during call ?
        angle_zy_set;                       // is antAngle_zy set during call ?

    // set-up grid
    scale           = 1;
    gridCfg.period  = 16*scale;
#if BOUNDARY == 1
    gridCfg.d_absorb= (int)(3*gridCfg.period);
#elif BOUNDARY == 2
    gridCfg.d_absorb= 8;
#endif
    gridCfg.Nx  = (400+0*200)*scale;
    gridCfg.Ny  = (300+0*100)*scale;
    gridCfg.Nz  = (200+0*150)*scale;
    gridCfg.Nz_ref  = 2*gridCfg.d_absorb + (int)gridCfg.period;
    gridCfg.t_end   = (int)((100-50)*gridCfg.period);

    gridCfg.B0_profile  = 0;
    gridCfg.ne_profile  = 2;

    // arrays realized as variable-length array (VLA)
    // E- and B-wavefield
    double (*EB_WAVE)[gridCfg.Ny][gridCfg.Nz]           = calloc(gridCfg.Nx, sizeof *EB_WAVE);
    double (*EB_WAVE_ref)[gridCfg.Ny][gridCfg.Nz_ref]   = calloc(gridCfg.Nx, sizeof *EB_WAVE_ref);
    // J-wavefield (in plasma) and background magnetic field
    double (*J_B0)[gridCfg.Ny][gridCfg.Nz]              = calloc(gridCfg.Nx, sizeof *J_B0);
    // background electron plasma density
    double (*n_e)[gridCfg.Ny/2][gridCfg.Nz/2]           = calloc(gridCfg.Nx/2, sizeof *n_e);
    // used when writing data into hdf5-files
    double (*data2save)[gridCfg.Ny/2][gridCfg.Nz/2]     = calloc(gridCfg.Nx/2, sizeof *data2save);
    // antenna: envelope of injected field
    double (*antField_xy)[gridCfg.Ny/2]                 = calloc(gridCfg.Nx/2, sizeof *antField_xy);
    // antenna: phase terms 
    double (*antPhaseTerms)[gridCfg.Ny/2]               = calloc(gridCfg.Nx/2, sizeof *antPhaseTerms);
    // time traces
    double (*timetraces)[8]                             = calloc((gridCfg.t_end/(int)gridCfg.period), sizeof *timetraces);

    // old E-fields required for Mur's boundary condition
#if BOUNDARY == 2
    double (*E_Xdir_OLD)[gridCfg.Ny][gridCfg.Nz]        = calloc(8,  sizeof *E_Xdir_OLD);
    double (*E_Ydir_OLD)[8][gridCfg.Nz]                 = calloc(gridCfg.Nx, sizeof *E_Ydir_OLD);
    double (*E_Zdir_OLD)[gridCfg.Ny][8]                 = calloc(gridCfg.Nx, sizeof *E_Zdir_OLD);
    double (*E_Xdir_OLD_ref)[gridCfg.Ny][gridCfg.Nz_ref]= calloc(8,  sizeof *E_Xdir_OLD_ref);
    double (*E_Ydir_OLD_ref)[8][gridCfg.Nz_ref]         = calloc(gridCfg.Nx, sizeof *E_Ydir_OLD_ref);
    double (*E_Zdir_OLD_ref)[gridCfg.Ny][8]             = calloc(gridCfg.Nx, sizeof *E_Zdir_OLD_ref);
#endif

    // array for detector antennas
    // sum_t(Ex*Ex) | sum_t(Ey*Ey) | sum_t(Ez*Ez) | sum_t(E*E) | rms(E)
#ifdef DETECTOR_ANTENNA_1D
    // TODO: change into 3D array, such that each detector antenna corresponds
    //       to one 2D array; that way it can be written much more failsafe...
    //       requires some changes in procedures for storing and saving
    double (*detAnt_01_fields)[5]       = calloc(gridCfg.Nx, sizeof *detAnt_01_fields);
    double (*detAnt_02_fields)[5]       = calloc(gridCfg.Nx, sizeof *detAnt_02_fields);
    double (*detAnt_03_fields)[5]       = calloc(gridCfg.Nx, sizeof *detAnt_03_fields);
    double (*detAnt_04_fields)[5]       = calloc(gridCfg.Nx, sizeof *detAnt_04_fields);
#endif

    // reading input parameter
    // used for checking if input parameter was provided
    angle_zx_set    = false;
    angle_zy_set    = false;
    // default values to be used if input parameter are not set
    beamCfg.antAngle_zx     = 0;
    beamCfg.antAngle_zy     = 0;
    // loop through input parameter
    printf( "number of input parameters provided during call: %d\n", argc-1 );
    while ( (opt_ret = getopt(argc, argv, "a:b:")) != -1 ){
        switch (opt_ret) {
            // angle between z=const plane and x=const plane
            case 'a': beamCfg.antAngle_zx   = atof(optarg);
                      angle_zx_set  = true;
                      break;
            case 'b': beamCfg.antAngle_zy   = atof(optarg);
                      angle_zy_set  = true;
                      break;
        }
    }
    if ( argc > 1 ) {
        printf( "following parameters were set during call: \n" );
        if (angle_zx_set)   printf( "    antAngle_zx = %f\n", beamCfg.antAngle_zx );
        if (angle_zy_set)   printf( "    antAngle_zy = %f\n", beamCfg.antAngle_zy );
    }

    beamCfg.exc_signal  = 5;//3;//4;
    beamCfg.rampUpMethod= 1;
    beamCfg.ant_x       = gridCfg.d_absorb + 8*gridCfg.period;//gridCfg.Nx/2;
    beamCfg.ant_y       = gridCfg.Ny/2;
    beamCfg.ant_z       = gridCfg.d_absorb + 4;
    // positions have to be even numbers, to ensure fields are accessed correctly
    if ((beamCfg.ant_x % 2) != 0)  ++beamCfg.ant_x;
    if ((beamCfg.ant_y % 2) != 0)  ++beamCfg.ant_y;
    if ((beamCfg.ant_z % 2) != 0)  ++beamCfg.ant_z;
    beamCfg.ant_w0x     = 2;
    beamCfg.ant_w0y     = 2;
    beamCfg.z2waist     = -(298.87)*.0;                // .2/l_0*period = -298.87

    pwr_dect    = gridCfg.d_absorb;

#ifdef DETECTOR_ANTENNA_1D
    detAnt_01_ypos  = beamCfg.ant_y;
    detAnt_01_zpos  = beamCfg.ant_z+2;
    detAnt_02_zpos  = round(beamCfg.ant_z+2 + 1*5*gridCfg.period); // steps of 5 cm for 28 GHz = 4.67*period
    detAnt_03_zpos  = round(beamCfg.ant_z+2 + 2*5*gridCfg.period);
    detAnt_04_zpos  = round(beamCfg.ant_z+2 + 3*5*gridCfg.period);
    // positions have to be even numbers, to ensure fields are accessed correctly
    if ((detAnt_01_ypos % 2) != 0)  ++detAnt_01_ypos;
    if ((detAnt_01_zpos % 2) != 0)  ++detAnt_01_zpos;
    if ((detAnt_02_zpos % 2) != 0)  ++detAnt_02_zpos;
    if ((detAnt_03_zpos % 2) != 0)  ++detAnt_03_zpos;
    if ((detAnt_04_zpos % 2) != 0)  ++detAnt_04_zpos;
    // issue a warning when detector antenna position is beyond Nz
    if (detAnt_04_zpos > (gridCfg.Nz - gridCfg.d_absorb)) {
        printf( "ERROR: check the detector antenna positions into z direction\n" );
        printf( "       Nz-d_absorb = %d, detAnt_04_zpos = %d\n", 
                gridCfg.Nz-gridCfg.d_absorb, detAnt_04_zpos );
    }
#endif

    // dt/dx = 0.5 is commenly used in 2D FDTD codes
    // Note that period refers to the wavelength in the numerical grid and not
    // in the "physical" grid (where one grid cell is equal to one Yee cell).
    // This means that in the physical grid, the wavelength is period/2, thus
    // in the equations we have to use period/2 for the wavelength.
    gridCfg.dx  = 1./(gridCfg.period/2);
    gridCfg.dt  = 1./(2.*(gridCfg.period/2));
        
#if BOUNDARY == 1
    eco         = 10./(double)(gridCfg.period);
#endif

    T_wave      = 0;
    omega_t     = .0;

    // the arrays are initialized with calloc() and thus don't require zeroing
    printf( "starting to set all variables to 0...\n" );
    power_abs_x1    = .0;
    power_abs_x2    = .0;
    power_abs_y1    = .0;
    power_abs_y2    = .0;
    power_abs_z1    = .0;
    power_abs_z2    = .0;
    power_abs_ref   = 1e-7;
    poynt_x1       = .0;
    poynt_x2       = .0;
    poynt_y1       = .0;
    poynt_y2       = .0;
    poynt_z1       = .0;
    poynt_z1_ref   = .0;
    poynt_z2       = .0;
    /*
    power_EE_x1    = .0;
    power_EE_x2    = .0;
    power_EE_y1    = .0;
    power_EE_y2    = .0;
    power_EE_z1    = .0;
    power_EE_z2    = .0;
    power_EE_ref   = .0;
    */
    printf( "...done setting all variables to 0\n" );

    printf( "starting do define antenna field...\n" );
    make_antenna_profile( &gridCfg, &beamCfg, 
                          antField_xy, antPhaseTerms );
    printf( "...done defining antenna field\n" );

    printf( "starting defining background plasma density\n" );
            // ne_profile: 1 = plasma mirror
            //             2 = linearly increasing profile
    make_density_profile( &gridCfg,  
            // cntrl_para: ne_profile=1 --> 0: plane mirror; oblique mirror: -.36397; 20 degrees: -.17633
            //             ne_profile=2 --> k0*Ln: 25
            25,
            n_e );
    printf( " ...setting density in absorber to 0...\n ");
    //set_densityInAbsorber_v2( &gridCfg, "z1", n_e );
    //set_densityInAbsorber_v2( &gridCfg, "x1x2y1y2z1", n_e );
    printf( "...done defining background plasma density\n" );

    printf( "starting defining background magnetic field...\n" );
    // B0x: even-odd-odd
    // B0y: odd-even-odd
    // B0z: odd-odd-even
            // B0_profile: 1 = constant field
    make_B0_profile(
            &gridCfg,
            // cntrl_para: B0_profile=1 --> value of Y
            .85, 
            J_B0 );
    printf( "...done defining background magnetic field\n" );

    // print some info to console
    printf( "Nx = %d, Ny = %d, Nz = %d\n", gridCfg.Nx, gridCfg.Ny, gridCfg.Nz );
    printf( "period = %d\n", (int)(gridCfg.period) );
    printf( "d_absorb = %d\n", gridCfg.d_absorb );
    printf( "t_end = %d\n", (int)(gridCfg.t_end) );
    printf( "antAngle_zx = %.2f, antAngle_zy = %.2f\n", beamCfg.antAngle_zx, beamCfg.antAngle_zy );
    printf( "ant_w0x = %.2f, ant_w0y = %.2f\n", beamCfg.ant_w0x, beamCfg.ant_w0y ); 
    printf( "ant_x = %d, ant_y = %d, ant_z = %d\n", beamCfg.ant_x, beamCfg.ant_y, beamCfg.ant_z );
    printf( "Boundary condition set to '%d'\n", BOUNDARY );
#ifdef DETECTOR_ANTENNA_1D
    printf( "detector antenna positions: z1 = %d, y1 = %d\n", detAnt_01_zpos, detAnt_01_ypos );
    printf( "detector antenna positions: z2 = %d, y1 = %d\n", detAnt_02_zpos, detAnt_01_ypos );
    printf( "detector antenna positions: z3 = %d, y1 = %d\n", detAnt_03_zpos, detAnt_01_ypos );
    printf( "detector antenna positions: z4 = %d, y1 = %d\n", detAnt_04_zpos, detAnt_01_ypos );
#endif

#ifdef _OPENMP
#pragma omp parallel private(n_threads)
    {
    n_threads = omp_get_num_threads();
    printf( "number of threads that will be used (OpenMP) = %d\n", n_threads );
    }
#endif


    for (t_int=0 ; t_int <=gridCfg.t_end ; ++t_int) {
        
        omega_t += 2.*M_PI/gridCfg.period;

        // to avoid precision problems when a lot of pi's are summed up        
        if (omega_t >= 2.*M_PI) {
            omega_t    += -2.*M_PI;
            T_wave     += 1;
            //printf("status: number of oscillation periods: %d (t_int= %d) \n",T_wave,t_int);
        }

        // add source
        add_source( &gridCfg, &beamCfg,
                    .85*0.,     // .85=Y, this values should be calculated/extracted from ne-profile
                    t_int, omega_t, 
                    antField_xy, antPhaseTerms, EB_WAVE );
        add_source_ref( &gridCfg, &beamCfg,
                        .85*0.,     // .85=Y, this values should be calculated/extracted from ne-profile
                        t_int, omega_t, 
                        antField_xy, antPhaseTerms, EB_WAVE_ref );

        // apply absorbers
#if BOUNDARY == 1
        apply_absorber(     &gridCfg, eco, EB_WAVE );
        apply_absorber_ref( &gridCfg, eco, EB_WAVE_ref );
#endif

        // advance J
        // Jx: odd-even-even
        // Jy: even-odd-even
        // Jz: even-even-odd
        // B0x: even-odd-odd
        // B0y: odd-even-odd
        // B0z: odd-odd-even
        advance_J( &gridCfg, EB_WAVE, J_B0, n_e );

        // advance B
        advance_B(     &gridCfg, EB_WAVE );
        advance_B_ref( &gridCfg, EB_WAVE_ref );
        
        // advance E
        advance_E(     &gridCfg, EB_WAVE,     J_B0 );
        advance_E_ref( &gridCfg, EB_WAVE_ref       );

        // optionally, apply numerical viscosity
        //apply_numerical_viscosity( &gridCfg, EB_WAVE );

        // apply Mur's boundary conditions
#if BOUNDARY == 2
        abc_Mur_1st( &gridCfg, "x1x2y1y2z1z2",  
                     EB_WAVE, E_Xdir_OLD, E_Ydir_OLD, E_Zdir_OLD );
        abc_Mur_1st_ref( &gridCfg, 
                         EB_WAVE_ref, E_Xdir_OLD_ref, E_Ydir_OLD_ref, E_Zdir_OLD_ref );
        abc_Mur_saveOldE_xdir(    &gridCfg, EB_WAVE, E_Xdir_OLD );
        abc_Mur_saveOldE_ydir(    &gridCfg, EB_WAVE, E_Ydir_OLD );
        abc_Mur_saveOldE_zdir(    &gridCfg, EB_WAVE, E_Zdir_OLD );
        abc_Mur_saveOldEref_xdir( &gridCfg, EB_WAVE_ref, E_Xdir_OLD_ref );
        abc_Mur_saveOldEref_ydir( &gridCfg, EB_WAVE_ref, E_Ydir_OLD_ref );
        abc_Mur_saveOldEref_zdir( &gridCfg, EB_WAVE_ref, E_Zdir_OLD_ref );
#endif

#ifdef DETECTOR_ANTENNA_1D
        // store wavefields for detector antennas over the final 10 
        // oscillation periods, it was found previously that only one period
        // does not result in a too nice average
        if ( t_int >= (gridCfg.t_end-10*gridCfg.period) ) {
            if (detAnt_01_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
                detAnt1D_storeValues( &gridCfg, detAnt_01_ypos, detAnt_01_zpos,
                                      t_int,  
                                      EB_WAVE, detAnt_01_fields );
            }
            if (detAnt_02_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
                detAnt1D_storeValues( &gridCfg, detAnt_01_ypos, detAnt_02_zpos,
                                      t_int, 
                                      EB_WAVE, detAnt_02_fields );
            }
            if (detAnt_03_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
                detAnt1D_storeValues( &gridCfg, detAnt_01_ypos, detAnt_03_zpos,
                                      t_int,
                                      EB_WAVE, detAnt_03_fields );
            }
            if (detAnt_04_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
                detAnt1D_storeValues( &gridCfg, detAnt_01_ypos, detAnt_04_zpos,
                                      t_int,
                                      EB_WAVE, detAnt_04_fields );
            }
        }
#endif

        // IQ detector for power detection
        if ( t_int >= 20*gridCfg.period ) {
            // z1-plane and z2-plane
            poynt_z1_ref    = calc_poynt_4( &gridCfg, pwr_dect, "ref_z1", EB_WAVE, EB_WAVE_ref );
            poynt_z1        = calc_poynt_4( &gridCfg, pwr_dect, "z1",     EB_WAVE, EB_WAVE_ref );
            poynt_z2        = calc_poynt_4( &gridCfg, pwr_dect, "z2",     EB_WAVE, EB_WAVE_ref );
            // x1-plane and x2-plane
            poynt_x1        = calc_poynt_4( &gridCfg, pwr_dect, "x1", EB_WAVE, EB_WAVE_ref );
            poynt_x2        = calc_poynt_4( &gridCfg, pwr_dect, "x2", EB_WAVE, EB_WAVE_ref );
            // y1-plane and y2-plane
            poynt_y1        = calc_poynt_4( &gridCfg, pwr_dect, "y1", EB_WAVE, EB_WAVE_ref );
            poynt_y2        = calc_poynt_4( &gridCfg, pwr_dect, "y2", EB_WAVE, EB_WAVE_ref );

            
//            printf( "t = %d, power_abs_ref = %13.5e, power_abs_z1 = %13.5e, power_abs_z2 = %13.5e, poynt_z1 = %13.5e, poynt_z2 = %13.5e\n",
//                    t_int, power_abs_ref, power_abs_z1, power_abs_z2, poynt_z1, poynt_z2 );

            power_abs_ref   = .99*power_abs_ref + .01*poynt_z1_ref;
            power_abs_z1    = .99*power_abs_z1  + .01*poynt_z1;
            power_abs_z2    = .99*power_abs_z2  + .01*poynt_z2;
            power_abs_x1    = .99*power_abs_x1  + .01*poynt_x1;
            power_abs_x2    = .99*power_abs_x2  + .01*poynt_x2;
            power_abs_y1    = .99*power_abs_y1  + .01*poynt_y1;
            power_abs_y2    = .99*power_abs_y2  + .01*poynt_y2;

            /*
            // EE
            // z1-plane and z2-plane
            power_EE_ref    += calc_power_EE_1( gridCfg.Nx, gridCfg.Ny, gridCfg.Nz, gridCfg.Nz_ref, gridCfg.d_absorb, "ref_z1", EB_WAVE, EB_WAVE_ref );
            power_EE_z1     += calc_power_EE_1( gridCfg.Nx, gridCfg.Ny, gridCfg.Nz, gridCfg.Nz_ref, gridCfg.d_absorb, "z1",     EB_WAVE, EB_WAVE_ref );
            power_EE_z2     += calc_power_EE_1( gridCfg.Nx, gridCfg.Ny, gridCfg.Nz, gridCfg.Nz_ref, gridCfg.d_absorb, "z2",     EB_WAVE, EB_WAVE_ref );
            // x1-plane and x2-plane
            power_EE_x1     += calc_power_EE_1( gridCfg.Nx, gridCfg.Ny, gridCfg.Nz, gridCfg.Nz_ref, gridCfg.d_absorb, "x1",     EB_WAVE, EB_WAVE_ref );
            power_EE_x2     += calc_power_EE_1( gridCfg.Nx, gridCfg.Ny, gridCfg.Nz, gridCfg.Nz_ref, gridCfg.d_absorb, "x2",     EB_WAVE, EB_WAVE_ref );
            // y1-plane and y2-plane
            power_EE_y1     += calc_power_EE_1( gridCfg.Nx, gridCfg.Ny, gridCfg.Nz, gridCfg.Nz_ref, gridCfg.d_absorb, "y1",     EB_WAVE, EB_WAVE_ref );
            power_EE_y2     += calc_power_EE_1( gridCfg.Nx, gridCfg.Ny, gridCfg.Nz, gridCfg.Nz_ref, gridCfg.d_absorb, "y2",     EB_WAVE, EB_WAVE_ref );
            */

        }


        if ( (t_int % (int)(gridCfg.period)) == 4 )  {
            printf( "status: number of oscillation periods: %d (t_int= %d) \n",T_wave,t_int);
            printf( "        Poynting-power: z1 = %13.6e, z2 = %13.6e, x1 = %13.6e, x2 = %13.6e, y1 = %13.6e, y2 = %13.6e, (z1+z2+x1+x2+y1+y2)/z1_ref = %13.6e %%\n",
                    power_abs_z1/power_abs_ref, 
                    power_abs_z2/power_abs_ref,
                    power_abs_x1/power_abs_ref, 
                    power_abs_x2/power_abs_ref,
                    power_abs_y1/power_abs_ref, 
                    power_abs_y2/power_abs_ref,
                    (power_abs_x1+power_abs_x2 + power_abs_y1+power_abs_y2 + power_abs_z1+power_abs_z2)/power_abs_ref * 100.
                    );
            /*
            printf( "        Power_EE_d-abs: z1 = %13.6e, z2 = %13.6e, x1 = %13.6e, x2 = %13.6e, y1 = %13.6e, y2 = %13.6e, ref = %13.6e\n",
                    power_EE_z1, 
                    power_EE_z2,
                    power_EE_x1, 
                    power_EE_x2,
                    power_EE_y1, 
                    power_EE_y2,
                    power_EE_ref
//                    (power_abs_x1+power_abs_x2 + power_abs_y1+power_abs_y2 + power_abs_z1+power_abs_z2)/power_abs_ref * 100.
                    );
            */
            timetraces[T_wave][0]   = (double)t_int;
            timetraces[T_wave][1]   = (double)T_wave;
            timetraces[T_wave][2]   = power_abs_z1/power_abs_ref;
            timetraces[T_wave][3]   = power_abs_z2/power_abs_ref;
            timetraces[T_wave][4]   = power_abs_x1/power_abs_ref;
            timetraces[T_wave][5]   = power_abs_x2/power_abs_ref;
            timetraces[T_wave][6]   = power_abs_y1/power_abs_ref;
            timetraces[T_wave][7]   = power_abs_y2/power_abs_ref;

        }
    } // end of time loop

    printf( "-------------------------------------------------------------------------------------------------------------\n" );
    printf( "  T   |   poynt_z1   |   poynt_z2   |   poynt_x1   |   poynt_x2   |   poynt_y1   |   poynt_y2   |  P_out     \n" );
    printf( "------+--------------+--------------+--------------+--------------+--------------+--------------+------------\n" );
    for ( ii=0 ; ii<(gridCfg.t_end/(int)gridCfg.period) ; ++ii )
        printf( " %4d |%13.6e |%13.6e |%13.6e |%13.6e |%13.6e |%13.6e |%13.6e\n",
                (int)timetraces[ii][1], //timetraces[ii][1],
                timetraces[ii][2], timetraces[ii][3],
                timetraces[ii][4], timetraces[ii][5],
                timetraces[ii][6], timetraces[ii][7],
                (timetraces[ii][2]+timetraces[ii][3] + timetraces[ii][4]+timetraces[ii][5] + timetraces[ii][6]+timetraces[ii][7])
              );
    printf( "-------------------------------------------------------------------------------------------------------------\n" );

    // write timetrace data into file
    // open file in w(rite) mode; might consider using a+ instead
    writeTimetraces2ascii( (gridCfg.t_end/(int)gridCfg.period), 8, gridCfg.t_end, gridCfg.period, 
                           "timetraces2.dat", timetraces );

    // save into hdf5
    // abs(E)
    // prepare array for that
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii<gridCfg.Nx ; ii+=2) {
        for (jj=0 ; jj<gridCfg.Ny ; jj+=2) {
            for (kk=0 ; kk<gridCfg.Nz ; kk+=2) {
                data2save[(ii/2)][(jj/2)][(kk/2)] = 
                    sqrt( pow(EB_WAVE[ii+1][jj  ][kk  ],2) 
                         +pow(EB_WAVE[ii  ][jj+1][kk  ],2) 
                         +pow(EB_WAVE[ii  ][jj  ][kk+1],2) );
            }
        }
    }
    len_str = snprintf( filename_hdf5, sizeof(filename_hdf5), "fileout.h5");
    if ( (len_str < 0) || (len_str >= sizeof(filename_hdf5)) ) {
        printf( "ERROR: could not write filename_hdf5 string\n" );  // use a proper error handler here
    } else {
        sprintf( dSet_name, "E_abs__tint%05d", t_int );
        printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( gridCfg.Nx/2, gridCfg.Ny/2, gridCfg.Nz/2, filename_hdf5, dSet_name, data2save) ) ;
    }
    set2zero_3D( gridCfg.Nx/2, gridCfg.Ny/2, gridCfg.Nz/2, data2save );
    // density
    sprintf( dSet_name, "n_e" );
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( gridCfg.Nx/2, gridCfg.Ny/2, gridCfg.Nz/2, filename_hdf5, dSet_name, n_e) ) ;
    // background magnetic field
    // B0x: even-odd-odd
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii<gridCfg.Nx ; ii+=2) {
        for (jj=0 ; jj<gridCfg.Ny ; jj+=2) {
            for (kk=0 ; kk<gridCfg.Nz ; kk+=2) {
                data2save[(ii/2)][(jj/2)][(kk/2)] = J_B0[ii  ][jj+1][kk+1];
            }
        }
    }
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( gridCfg.Nx/2, gridCfg.Ny/2, gridCfg.Nz/2, filename_hdf5, "B0x", data2save) ) ;
    set2zero_3D( gridCfg.Nx/2, gridCfg.Ny/2, gridCfg.Nz/2, data2save );
    // B0y: odd-even-odd
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii<gridCfg.Nx ; ii+=2) {
        for (jj=0 ; jj<gridCfg.Ny ; jj+=2) {
            for (kk=0 ; kk<gridCfg.Nz ; kk+=2) {
                data2save[(ii/2)][(jj/2)][(kk/2)] = J_B0[ii+1][jj  ][kk+1];
            }
        }
    }
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( gridCfg.Nx/2, gridCfg.Ny/2, gridCfg.Nz/2, filename_hdf5, "B0y", data2save) ) ;
    set2zero_3D( gridCfg.Nx/2, gridCfg.Ny/2, gridCfg.Nz/2, data2save );
    // B0z: odd-odd-even
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii<gridCfg.Nx ; ii+=2) {
        for (jj=0 ; jj<gridCfg.Ny ; jj+=2) {
            for (kk=0 ; kk<gridCfg.Nz ; kk+=2) {
                data2save[(ii/2)][(jj/2)][(kk/2)] = J_B0[ii+1][jj+1][kk  ];
            }
        }
    }
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( gridCfg.Nx/2, gridCfg.Ny/2, gridCfg.Nz/2, filename_hdf5, "B0z", data2save) ) ;
    set2zero_3D( gridCfg.Nx/2, gridCfg.Ny/2, gridCfg.Nz/2, data2save );

    writeConfig2HDF( &gridCfg, &beamCfg, filename_hdf5 );


#if defined(HDF5) && defined(DETECTOR_ANTENNA_1D)
    if (detAnt_01_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
        detAnt1D_write2hdf5( gridCfg.Nx, filename_hdf5, "/detAnt_01" , 
                             detAnt_01_ypos, detAnt_01_zpos,
                             detAnt_01_fields );
    }
    if (detAnt_02_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
        detAnt1D_write2hdf5( gridCfg.Nx, filename_hdf5, "/detAnt_02" , 
                             detAnt_01_ypos, detAnt_02_zpos,
                             detAnt_02_fields );
    }
    if (detAnt_03_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
        detAnt1D_write2hdf5( gridCfg.Nx, filename_hdf5, "/detAnt_03" , 
                             detAnt_01_ypos, detAnt_03_zpos,
                             detAnt_03_fields );
    }
    if (detAnt_04_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
        detAnt1D_write2hdf5( gridCfg.Nx, filename_hdf5, "/detAnt_04" , 
                             detAnt_01_ypos, detAnt_04_zpos,
                             detAnt_04_fields );
    }
#endif

    free( EB_WAVE );
    printf( "freed EB_WAVE\n" );
    free( J_B0 );
    printf( "freed J_B0\n" );
    free( n_e );
    printf( "freed n_e\n" );
    free( data2save );
    printf( "freed data2save\n" );
    return EXIT_SUCCESS;
}//}}}


int make_density_profile( gridConfiguration *gridCfg, 
                          double cntrl_para, 
                          double n_e[gridCfg->Nx/2][gridCfg->Ny/2][gridCfg->Nz/2] ) {
//{{{
    size_t
        ii, jj, kk, 
        ne_start_z;
    double
        ne_max,
        ne_k0Ln,
        aux;

    // if density is larger than this value, FDTD code becomes instable
    ne_max  = gridCfg->period * 2./5.;

    if ( gridCfg->ne_profile == 1 ) {
        // plasma mirror
        for (ii=0 ; ii<(gridCfg->Nx/2) ; ++ii) {
            for (jj=0 ; jj<(gridCfg->Ny/2) ; ++jj) {
                //for (kk=((gridCfg->Nz-gridCfg->d_absorb-4)/2) ; kk<(gridCfg->Nz/2) ; ++kk) {
                for (kk=0 ; kk<(gridCfg->Nz/2) ; ++kk) {
                    // z = m*y + b = delta_z/delta_y * y + z_start
                    //             = cntrl_para * y + z_start
                    if ( (double)kk > (cntrl_para*(double)jj + ((double)gridCfg->Nz-(double)gridCfg->d_absorb-4.)/2.) ) {
                        n_e[ii][jj][kk] = ne_max;
                    }
                }
            }
        }
    } else if ( gridCfg->ne_profile == 2 ) {
        // linearly increasing profile with k0Ln as slope
        // n_e(z) = m*z 
        // with m = 2*pi/(k0Ln*lambda) 
        //      z = z-starting_position
        ne_start_z  = (gridCfg->d_absorb + gridCfg->period)/2;
        //ne_start_z  = (gridCfg->d_absorb + 224.155)/2;    // benchmark scenario from STEP project: .15m/l_0*period
        if (ne_start_z%2 != 0)
            ne_start_z  += 1;
        ne_k0Ln     = cntrl_para;
        printf( "make_density_profile: ne_profile = %d, ne_start_z = %ld, k0Ln = %f\n", 
                gridCfg->ne_profile, ne_start_z, ne_k0Ln );
        for (ii=0 ; ii<(gridCfg->Nx/2) ; ++ii) {
            for (jj=0 ; jj<(gridCfg->Ny/2) ; ++jj) {
                for (kk=0 ; kk<(gridCfg->Nz/2) ; ++kk) {
                    aux = ((double)kk - (double)ne_start_z) * (2.*M_PI / (ne_k0Ln*gridCfg->period));
                    // negative density values are unphysical
                    if (aux < .0)
                        aux = .0;
                    // FDTD full-wave codes only allow for a maximum density value
                    if (aux > ne_max) {
                        aux  = ne_max;
                        printf( "    maximum density achieved (ii, jj, kk = %ld, %ld, %ld): %f\n", ii, jj, kk, aux );
                    }
                    //if ((ii%10 == 0) && (jj%10 == 0) && (kk%10 == 0))
                    //if (kk%10 == 0)
                    //    printf( "    n_e[%4ld][%4ld][%4ld] = %f\n", ii, jj, kk, aux );
                    n_e[ii][jj][kk] = aux;
                }
            }
        }
    } else if ( gridCfg->ne_profile == 3 ) {
        // plasma cylinder (2D gauss) in center of grid
        //
        //                               (y-y0)^2        (z-z0)^2
        // n_e(y,z) = n_e,max * exp( -( ------------ + ------------- ) )
        //                               2*sigma_y^2    2*sigma_z^2
        //
        for (ii=0 ; ii<(gridCfg->Nx/2) ; ++ii) {
            for (jj=0 ; jj<(gridCfg->Ny/2) ; ++jj) {
                for (kk=0 ; kk<(gridCfg->Nz/2) ; ++kk) {
                    n_e[ii][jj][kk]    = exp( -1.* (
                                                 pow((double)jj-(double)gridCfg->Ny/4., 2)/(2*pow(gridCfg->period/2.,2)) 
                                                +pow((double)kk-(double)gridCfg->Nz/4., 2)/(2*pow(gridCfg->period/2.,2))
                                             )) * 5.;
                }
            }
        }
    } else if ( gridCfg->ne_profile == 4 ) {
        // same as ne_profile = 3, but plasma cylinder is now along y
        for (ii=0 ; ii<(gridCfg->Nx/2) ; ++ii) {
            for (jj=0 ; jj<(gridCfg->Ny/2) ; ++jj) {
                for (kk=0 ; kk<(gridCfg->Nz/2) ; ++kk) {
                    n_e[ii][jj][kk]    = exp( -1.* (
                                                 pow((double)ii-(double)gridCfg->Nx/4., 2)/(2*pow(gridCfg->period/2.,2)) 
                                                +pow((double)kk-(double)gridCfg->Nz/4., 2)/(2*pow(gridCfg->period/2.,2))
                                             )) * 2;//5.;
                }
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int set_densityInAbsorber_v2( gridConfiguration *gridCfg, 
                              char absorber[], 
                              double n_e[gridCfg->Nx/2][gridCfg->Ny/2][gridCfg->Nz/2] ) {
//{{{

    double
        x, y, z,
        x0,x1, y0,y1, z0,z1,
        ne_absorb,              // electron density in absorber
        smooth,                 // 1: relatively steep, .2: more smooth
        ne_dist,                // proportional to distance to absorber, where n_e starts to decrease
        scale_fact;

    ne_absorb   = .0;
    smooth      = .5;//.2;
    ne_dist     = round( gridCfg->period/1 );

    x0          = (double)gridCfg->d_absorb + ne_dist;
    x1          = (double)gridCfg->Nx - (gridCfg->d_absorb + ne_dist);
    y0          = (double)gridCfg->d_absorb + ne_dist;
    y1          = (double)gridCfg->Ny - (gridCfg->d_absorb + ne_dist);
    z0          = (double)gridCfg->d_absorb + ne_dist;
    z1          = (double)gridCfg->Nz - (gridCfg->d_absorb + ne_dist);

    // scale to density grid which is only half the size of FDTD-wavefields grid
    // since 2 numerical grid points correspond to one "physical" grid point
    // and we want not to vary the background parameters within one physical
    // grid point
    x0  = round(x0*.5);
    x1  = round(x1*.5);
    y0  = round(y0*.5);
    y1  = round(y1*.5);
    z0  = round(z0*.5);
    z1  = round(z1*.5);

    // the string "absorber" is used to set in which absorber n_e will be modified
    // the comparison is done with the strstr() function, which return the address
    // of the substring if found, NULL otherwise

    // set density in x0 absorber
    if ( strstr(absorber,"x1") ) {
        for ( x=0; x<(gridCfg->Nx/2) ; ++x ) {
            scale_fact  = +.5*(    tanh(smooth*(x-x0)) + 1);        // x0 boundary
            //printf( "x1: x=%.1f, scale_fact=%f\n", x, scale_fact) ;
            for ( y=0. ; y<(gridCfg->Ny/2) ; ++y )   {
                for (z=0 ; z<(gridCfg->Nz/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }
    // set density in x1 absorber
    if ( strstr(absorber,"x2") ) {
        for ( x=0; x<(gridCfg->Nx/2) ; ++x ) {
            scale_fact  = +.5*(-1.*tanh(smooth*(x-x1)) + 1);       // x1 boundary
            //printf( "x2: x=%.1f, scale_fact=%f\n", x, scale_fact) ;
            for ( y=0. ; y<(gridCfg->Ny/2) ; ++y )   {
                for (z=0 ; z<(gridCfg->Nz/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }

    // set density in y0 absorber
    if ( strstr(absorber,"y1") ) {
        for ( y=0; y<(gridCfg->Ny/2) ; ++y ) {
            scale_fact  = +.5*(    tanh(smooth*(y-y0)) + 1);        // y0 boundary
            //printf( "y1: y=%.1f, scale_fact=%f\n", y, scale_fact) ;
            for ( x=0; x<(gridCfg->Nx/2) ; ++x ) {
                for (z=0 ; z<(gridCfg->Nz/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }
    // set density in y1 absorber
    if ( strstr(absorber,"y2") ) {
        for ( y=0; y<(gridCfg->Ny/2) ; ++y ) {
            scale_fact  = +.5*(-1.*tanh(smooth*(y-y1)) + 1);       // y1 boundary
            //printf( "y2: y=%.1f, scale_fact=%f\n", y, scale_fact) ;
            for ( x=0; x<(gridCfg->Nx/2) ; ++x ) {
                for (z=0 ; z<(gridCfg->Nz/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }

    // set density in z0 absorber
    if ( strstr(absorber,"z1") ) {
        for ( z=0 ; z<(gridCfg->Nz/2) ; ++z) {
            scale_fact  = +.5*(    tanh(smooth*(z-z0)) + 1);        // z0 boundary
            //printf( "z1: z=%.1f, scale_fact=%f\n", z, scale_fact) ;
            for ( x=0; x<(gridCfg->Nx/2) ; ++x ) {
                for ( y=0; y<(gridCfg->Ny/2) ; ++y ) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }
    // set density in z1 absorber
    if ( strstr(absorber,"z2") ) {
        for ( z=0 ; z<(gridCfg->Nz/2) ; ++z) {
            scale_fact  = +.5*(-1.*tanh(smooth*(z-z1)) + 1);       // z1 boundary
            //printf( "z2: z=%.1f, scale_fact=%f\n", z, scale_fact) ;
            for ( x=0; x<(gridCfg->Nx/2) ; ++x ) {
                for ( y=0; y<(gridCfg->Ny/2) ; ++y ) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }

    return EXIT_SUCCESS;
} //}}}


int make_B0_profile( gridConfiguration *gridCfg, 
                     double cntrl_para, 
                     double J_B0[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] ) {
//{{{
    size_t
        ii, jj, kk; 
//    double
//        aux;

    // B0x: even-odd-odd
    // B0y: odd-even-odd
    // B0z: odd-odd-even
    if ( gridCfg->B0_profile == 1 ) {
        // constant field
        for (ii=0 ; ii<gridCfg->Nx ; ii+=2) {
            for (jj=0 ; jj<gridCfg->Ny ; jj+=2) {
                for (kk=0 ; kk<gridCfg->Nz ; kk+=2) {
                    J_B0[ii  ][jj+1][kk+1] = cntrl_para;
                    J_B0[ii+1][jj  ][kk+1] = cntrl_para*.0;
                    J_B0[ii+1][jj+1][kk  ] = cntrl_para*.0;
                }
            }
        }
    } 
    return EXIT_SUCCESS;
}//}}}


int apply_absorber( gridConfiguration *gridCfg, 
                    double eco, 
                    double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] ) {
//{{{
    size_t
        ii, jj, kk;
    double
        damp;

    // z1 absorber: z=0...d_absorb
//#pragma omp parallel for collapse(2) default(shared) private(k,j,damp) <-- can collapse be used here?
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->d_absorb-2 ; kk+=2) {
                damp = ((double)kk-(double)gridCfg->d_absorb)/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
//                if ((ii%10 == 0) && (jj%10 == 0) && (kk%10 == 0)) 
//                    printf( "z1: ii=%3d, jj=%3d, kk=%3d, (kk-d_abs)/d_abs=%f, damp=%f\n", 
//                            ii, jj, kk, ((double)kk-(double)d_absorb)/(double)d_absorb, damp );
            }
        }
    }
    // z2 absorber: z=d_absorb...Nz
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=(gridCfg->Nz-gridCfg->d_absorb) ; kk<gridCfg->Nz-2 ; kk+=2) {      //Nz-d_absorb-2 ???
                damp = ((double)kk-((double)gridCfg->Nz-(double)gridCfg->d_absorb))/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }      
    // x1 absorber: x=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
            for (ii=2 ; ii<gridCfg->d_absorb-2 ; ii+=2) {
                damp = ((double)ii-(double)gridCfg->d_absorb)/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // x2 absorber: x=d_absorb...Nx
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {  
            for (ii=(gridCfg->Nx-gridCfg->d_absorb) ; ii<gridCfg->Nx-2 ; ii+=2) {    //Nx-d_absorb-2 ???
                damp = ((double)ii-((double)gridCfg->Nx-(double)gridCfg->d_absorb))/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // y1 absorber: y=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
            for (jj=2 ; jj<gridCfg->d_absorb-2 ; jj+=2) {
                damp = ((double)jj-(double)gridCfg->d_absorb)/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // y2 absorber: y=d_absorb...Ny
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
            for (jj=(gridCfg->Ny-gridCfg->d_absorb) ; jj<gridCfg->Ny-2 ; jj+=2) {  //Ny-d_absorb-2 ???
                damp = ((double)jj-((double)gridCfg->Ny-(double)gridCfg->d_absorb))/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int apply_absorber_ref( gridConfiguration *gridCfg, 
                        double eco, 
                        double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref] ) {
//{{{
    size_t
        ii, jj, kk;
    double
        damp;

    // z1 absorber: z=0...d_absorb
//#pragma omp parallel for collapse(2) default(shared) private(k,j,damp) <-- can collapse be used here?
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->d_absorb-2 ; kk+=2) {
                damp = ((double)kk-(double)gridCfg->d_absorb)/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
//                if ((ii%10 == 0) && (jj%10 == 0) && (kk%10 == 0)) 
//                    printf( "z1: ii=%3d, jj=%3d, kk=%3d, (kk-d_abs)/d_abs=%f, damp=%f\n", 
//                            ii, jj, kk, ((double)kk-(double)d_absorb)/(double)d_absorb, damp );
            }
        }
    }
    // z2 absorber: z=d_absorb...Nz
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=(gridCfg->Nz_ref-gridCfg->d_absorb) ; kk<gridCfg->Nz_ref-2 ; kk+=2) {      //Nz-d_absorb-2 ???
                damp = ((double)kk-((double)gridCfg->Nz_ref-(double)gridCfg->d_absorb))/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }      
    // x1 absorber: x=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            for (ii=2 ; ii<gridCfg->d_absorb-2 ; ii+=2) {
                damp = ((double)ii-(double)gridCfg->d_absorb)/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // x2 absorber: x=d_absorb...Nx
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {  
            for (ii=(gridCfg->Nx-gridCfg->d_absorb) ; ii<gridCfg->Nx-2 ; ii+=2) {    //Nx-d_absorb-2 ???
                damp = ((double)ii-((double)gridCfg->Nx-(double)gridCfg->d_absorb))/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // y1 absorber: y=0...d_absorb
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            for (jj=2 ; jj<gridCfg->d_absorb-2 ; jj+=2) {
                damp = ((double)jj-(double)gridCfg->d_absorb)/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    // y2 absorber: y=d_absorb...Ny
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            for (jj=(gridCfg->Ny-gridCfg->d_absorb) ; jj<gridCfg->Ny-2 ; jj+=2) {  //Ny-d_absorb-2 ???
                damp = ((double)jj-((double)gridCfg->Ny-(double)gridCfg->d_absorb))/(double)gridCfg->d_absorb;
                damp = ABSORBER_DAMPING(eco,damp);

                EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                EB_WAVE[ii  ][jj  ][kk+1] *= damp;
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int apply_absorber_v2( size_t N_x, size_t N_y, size_t N_z, int d_absorb, double eco, 
                       char absorber[],
                       double EB_WAVE[N_x][N_y][N_z] ) {
//{{{
    size_t
        ii, jj, kk;
    double
        damp;

    // z1 absorber: z=0...d_absorb
//#pragma omp parallel for collapse(2) default(shared) private(k,j,damp) <-- can collapse be used here?
    if ( strstr(absorber,"z1") ) {      
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (jj=2 ; jj<N_y-2 ; jj+=2) {
                for (kk=2 ; kk<d_absorb-2 ; kk+=2) {
                    damp = ((double)kk-(double)d_absorb)/(double)d_absorb;
                    damp = ABSORBER_DAMPING(eco,damp);

                    EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                    EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                    EB_WAVE[ii  ][jj  ][kk+1] *= damp;
//                if ((ii%10 == 0) && (jj%10 == 0) && (kk%10 == 0)) 
//                    printf( "z1: ii=%3d, jj=%3d, kk=%3d, (kk-d_abs)/d_abs=%f, damp=%f\n", 
//                            ii, jj, kk, ((double)kk-(double)d_absorb)/(double)d_absorb, damp );
                }
            }
        }
    }
    // z2 absorber: z=d_absorb...Nz
    if ( strstr(absorber,"z2") ) {      
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (jj=2 ; jj<N_y-2 ; jj+=2) {
                for (kk=(N_z-d_absorb) ; kk<N_z-2 ; kk+=2) {      //Nz-d_absorb-2 ???
                    damp = ((double)kk-((double)N_z-(double)d_absorb))/(double)d_absorb;
                    damp = ABSORBER_DAMPING(eco,damp);

                    EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                    EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                    EB_WAVE[ii  ][jj  ][kk+1] *= damp;
                }
            }
        }      
    }
    // x1 absorber: x=0...d_absorb
    if ( strstr(absorber,"x1") ) {
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (jj=2 ; jj<N_y-2 ; jj+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {
                for (ii=2 ; ii<d_absorb-2 ; ii+=2) {
                    damp = ((double)ii-(double)d_absorb)/(double)d_absorb;
                    damp = ABSORBER_DAMPING(eco,damp);

                    EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                    EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                    EB_WAVE[ii  ][jj  ][kk+1] *= damp;
                }
            }
        }
    }
    if ( strstr(absorber,"x2") ) {
    // x2 absorber: x=d_absorb...Nx
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (jj=2 ; jj<N_y-2 ; jj+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {  
                for (ii=(N_x-d_absorb) ; ii<N_x-2 ; ii+=2) {    //Nx-d_absorb-2 ???
                    damp = ((double)ii-((double)N_x-(double)d_absorb))/(double)d_absorb;
                    damp = ABSORBER_DAMPING(eco,damp);

                    EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                    EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                    EB_WAVE[ii  ][jj  ][kk+1] *= damp;
                }
            }
        }
    }
    // y1 absorber: y=0...d_absorb
    if ( strstr(absorber,"y1") ) {
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {
                for (jj=2 ; jj<d_absorb-2 ; jj+=2) {
                    damp = ((double)jj-(double)d_absorb)/(double)d_absorb;
                    damp = ABSORBER_DAMPING(eco,damp);

                    EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                    EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                    EB_WAVE[ii  ][jj  ][kk+1] *= damp;
                }
            }
        }
    }
    // y2 absorber: y=d_absorb...Ny
    if ( strstr(absorber,"y2") ) {
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {
                for (jj=(N_y-d_absorb) ; jj<N_y-2 ; jj+=2) {  //Ny-d_absorb-2 ???
                    damp = ((double)jj-((double)N_y-(double)d_absorb))/(double)d_absorb;
                    damp = ABSORBER_DAMPING(eco,damp);

                    EB_WAVE[ii+1][jj  ][kk  ] *= damp;
                    EB_WAVE[ii  ][jj+1][kk  ] *= damp;
                    EB_WAVE[ii  ][jj  ][kk+1] *= damp;
                }
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int apply_numerical_viscosity( gridConfiguration *gridCfg,
                               double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] ) {
    //{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk;

    double
        aux, ny;

    ny  = 1e-4;
    
#pragma omp parallel for default(shared) private(ii,jj,kk,aux)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
                aux = ny*(   EB_WAVE[ii  ][jj  ][kk+1+2] 
                          +  EB_WAVE[ii  ][jj  ][kk+1-2]
                          -2*EB_WAVE[ii  ][jj  ][kk+1  ]
                        );
                EB_WAVE[ii  ][jj  ][kk+1] += aux;
            }
        }
    }

    return EXIT_SUCCESS;
}//}}}


int abc_Mur_saveOldE_xdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[8][gridCfg->Ny][gridCfg->Nz] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        jj, kk, 
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
            // store values at x=0 and x=1
            // Ex: odd-even-even
            E_old[0+1][jj  ][kk  ]  = EB_WAVE[0+offset+1][jj  ][kk  ];
            E_old[2+1][jj  ][kk  ]  = EB_WAVE[2+offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_old[0  ][jj+1][kk  ]  = EB_WAVE[0+offset  ][jj+1][kk  ];
            E_old[2  ][jj+1][kk  ]  = EB_WAVE[2+offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_old[0  ][jj  ][kk+1]  = EB_WAVE[0+offset  ][jj  ][kk+1];
            E_old[2  ][jj  ][kk+1]  = EB_WAVE[2+offset  ][jj  ][kk+1];

            // store values at x=Nx-1 and x=Nx-2
            // Ex: odd-even-even
            E_old[4+1][jj  ][kk  ]  = EB_WAVE[gridCfg->Nx-4-offset+1][jj  ][kk  ];
            E_old[6+1][jj  ][kk  ]  = EB_WAVE[gridCfg->Nx-2-offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_old[4  ][jj+1][kk  ]  = EB_WAVE[gridCfg->Nx-4-offset  ][jj+1][kk  ];
            E_old[6  ][jj+1][kk  ]  = EB_WAVE[gridCfg->Nx-2-offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_old[4  ][jj  ][kk+1]  = EB_WAVE[gridCfg->Nx-4-offset  ][jj  ][kk+1];
            E_old[6  ][jj  ][kk+1]  = EB_WAVE[gridCfg->Nx-2-offset  ][jj  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldE_ydir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[gridCfg->Nx][8][gridCfg->Nz] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, kk,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
            // store values at y=0 and y=1
            // Ex: odd-even-even
            E_old[ii+1][0  ][kk  ]  = EB_WAVE[ii+1][0+offset  ][kk  ];
            E_old[ii+1][2  ][kk  ]  = EB_WAVE[ii+1][2+offset  ][kk  ];
            // Ey: even-odd-even
            E_old[ii  ][0+1][kk  ]  = EB_WAVE[ii  ][0+offset+1][kk  ];
            E_old[ii  ][2+1][kk  ]  = EB_WAVE[ii  ][2+offset+1][kk  ];
            // Ez: even-even-odd
            E_old[ii  ][0  ][kk+1]  = EB_WAVE[ii  ][0+offset  ][kk+1];
            E_old[ii  ][2  ][kk+1]  = EB_WAVE[ii  ][2+offset  ][kk+1];

            // store values at x=Nx-1 and x=Nx-2
            // Ex: odd-even-even
            E_old[ii+1][4  ][kk  ]  = EB_WAVE[ii+1][gridCfg->Ny-4-offset  ][kk  ];
            E_old[ii+1][6  ][kk  ]  = EB_WAVE[ii+1][gridCfg->Ny-2-offset  ][kk  ];
            // Ey: even-odd-even
            E_old[ii  ][4+1][kk  ]  = EB_WAVE[ii  ][gridCfg->Ny-4-offset+1][kk  ];
            E_old[ii  ][6+1][kk  ]  = EB_WAVE[ii  ][gridCfg->Ny-2-offset+1][kk  ];
            // Ez: even-even-odd
            E_old[ii  ][4  ][kk+1]  = EB_WAVE[ii  ][gridCfg->Ny-4-offset  ][kk+1];
            E_old[ii  ][6  ][kk+1]  = EB_WAVE[ii  ][gridCfg->Ny-2-offset  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldE_zdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[gridCfg->Nx][gridCfg->Ny][8] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            // store values at z=0 and z=1
            // Ex: odd-even-even
            E_old[ii+1][jj  ][0  ]  = EB_WAVE[ii+1][jj  ][0+offset  ];
            E_old[ii+1][jj  ][2  ]  = EB_WAVE[ii+1][jj  ][2+offset  ];
            // Ey: even-odd-even
            E_old[ii  ][jj+1][0  ]  = EB_WAVE[ii  ][jj+1][0+offset  ];
            E_old[ii  ][jj+1][2  ]  = EB_WAVE[ii  ][jj+1][2+offset  ];
            // Ez: even-even-odd
            E_old[ii  ][jj  ][0+1]  = EB_WAVE[ii  ][jj  ][0+offset+1];
            E_old[ii  ][jj  ][2+1]  = EB_WAVE[ii  ][jj  ][2+offset+1];

            // store values at z=Nz-1 and z=Nz-2
            // Ex: odd-even-even
            E_old[ii+1][jj  ][4  ]  = EB_WAVE[ii+1][jj  ][gridCfg->Nz-4-offset  ];
            E_old[ii+1][jj  ][6  ]  = EB_WAVE[ii+1][jj  ][gridCfg->Nz-2-offset  ];
            // Ey: even-odd-even
            E_old[ii  ][jj+1][4  ]  = EB_WAVE[ii  ][jj+1][gridCfg->Nz-4-offset  ];
            E_old[ii  ][jj+1][6  ]  = EB_WAVE[ii  ][jj+1][gridCfg->Nz-2-offset  ];
            // Ez: even-even-odd
            E_old[ii  ][jj  ][4+1]  = EB_WAVE[ii  ][jj  ][gridCfg->Nz-4-offset+1];
            E_old[ii  ][jj  ][6+1]  = EB_WAVE[ii  ][jj  ][gridCfg->Nz-2-offset+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldEref_xdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[8][gridCfg->Ny][gridCfg->Nz_ref] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        jj, kk, 
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            // store values at x=0 and x=1
            // Ex: odd-even-even
            E_old[0+1][jj  ][kk  ]  = EB_WAVE_ref[0+offset+1][jj  ][kk  ];
            E_old[2+1][jj  ][kk  ]  = EB_WAVE_ref[2+offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_old[0  ][jj+1][kk  ]  = EB_WAVE_ref[0+offset  ][jj+1][kk  ];
            E_old[2  ][jj+1][kk  ]  = EB_WAVE_ref[2+offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_old[0  ][jj  ][kk+1]  = EB_WAVE_ref[0+offset  ][jj  ][kk+1];
            E_old[2  ][jj  ][kk+1]  = EB_WAVE_ref[2+offset  ][jj  ][kk+1];

            // store values at x=Nx-1 and x=Nx-2
            // Ex: odd-even-even
            E_old[4+1][jj  ][kk  ]  = EB_WAVE_ref[gridCfg->Nx-4-offset+1][jj  ][kk  ];
            E_old[6+1][jj  ][kk  ]  = EB_WAVE_ref[gridCfg->Nx-2-offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_old[4  ][jj+1][kk  ]  = EB_WAVE_ref[gridCfg->Nx-4-offset  ][jj+1][kk  ];
            E_old[6  ][jj+1][kk  ]  = EB_WAVE_ref[gridCfg->Nx-2-offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_old[4  ][jj  ][kk+1]  = EB_WAVE_ref[gridCfg->Nx-4-offset  ][jj  ][kk+1];
            E_old[6  ][jj  ][kk+1]  = EB_WAVE_ref[gridCfg->Nx-2-offset  ][jj  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldEref_ydir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[gridCfg->Nx][8][gridCfg->Nz_ref] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, kk,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            // store values at y=0 and y=1
            // Ex: odd-even-even
            E_old[ii+1][0  ][kk  ]  = EB_WAVE_ref[ii+1][0+offset  ][kk  ];
            E_old[ii+1][2  ][kk  ]  = EB_WAVE_ref[ii+1][2+offset  ][kk  ];
            // Ey: even-odd-even
            E_old[ii  ][0+1][kk  ]  = EB_WAVE_ref[ii  ][0+offset+1][kk  ];
            E_old[ii  ][2+1][kk  ]  = EB_WAVE_ref[ii  ][2+offset+1][kk  ];
            // Ez: even-even-odd
            E_old[ii  ][0  ][kk+1]  = EB_WAVE_ref[ii  ][0+offset  ][kk+1];
            E_old[ii  ][2  ][kk+1]  = EB_WAVE_ref[ii  ][2+offset  ][kk+1];

            // store values at x=Nx-1 and x=Nx-2
            // Ex: odd-even-even
            E_old[ii+1][4  ][kk  ]  = EB_WAVE_ref[ii+1][gridCfg->Ny-4-offset  ][kk  ];
            E_old[ii+1][6  ][kk  ]  = EB_WAVE_ref[ii+1][gridCfg->Ny-2-offset  ][kk  ];
            // Ey: even-odd-even
            E_old[ii  ][4+1][kk  ]  = EB_WAVE_ref[ii  ][gridCfg->Ny-4-offset+1][kk  ];
            E_old[ii  ][6+1][kk  ]  = EB_WAVE_ref[ii  ][gridCfg->Ny-2-offset+1][kk  ];
            // Ez: even-even-odd
            E_old[ii  ][4  ][kk+1]  = EB_WAVE_ref[ii  ][gridCfg->Ny-4-offset  ][kk+1];
            E_old[ii  ][6  ][kk+1]  = EB_WAVE_ref[ii  ][gridCfg->Ny-2-offset  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldEref_zdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[gridCfg->Nx][gridCfg->Ny][8] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            // store values at z=0 and z=1
            // Ex: odd-even-even
            E_old[ii+1][jj  ][0  ]  = EB_WAVE_ref[ii+1][jj  ][0+offset  ];
            E_old[ii+1][jj  ][2  ]  = EB_WAVE_ref[ii+1][jj  ][2+offset  ];
            // Ey: even-odd-even
            E_old[ii  ][jj+1][0  ]  = EB_WAVE_ref[ii  ][jj+1][0+offset  ];
            E_old[ii  ][jj+1][2  ]  = EB_WAVE_ref[ii  ][jj+1][2+offset  ];
            // Ez: even-even-odd
            E_old[ii  ][jj  ][0+1]  = EB_WAVE_ref[ii  ][jj  ][0+offset+1];
            E_old[ii  ][jj  ][2+1]  = EB_WAVE_ref[ii  ][jj  ][2+offset+1];

            // store values at z=Nz-1 and z=Nz-2
            // Ex: odd-even-even
            E_old[ii+1][jj  ][4  ]  = EB_WAVE_ref[ii+1][jj  ][gridCfg->Nz_ref-4-offset  ];
            E_old[ii+1][jj  ][6  ]  = EB_WAVE_ref[ii+1][jj  ][gridCfg->Nz_ref-2-offset  ];
            // Ey: even-odd-even
            E_old[ii  ][jj+1][4  ]  = EB_WAVE_ref[ii  ][jj+1][gridCfg->Nz_ref-4-offset  ];
            E_old[ii  ][jj+1][6  ]  = EB_WAVE_ref[ii  ][jj+1][gridCfg->Nz_ref-2-offset  ];
            // Ez: even-even-odd
            E_old[ii  ][jj  ][4+1]  = EB_WAVE_ref[ii  ][jj  ][gridCfg->Nz_ref-4-offset+1];
            E_old[ii  ][jj  ][6+1]  = EB_WAVE_ref[ii  ][jj  ][gridCfg->Nz_ref-2-offset+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_1st( gridConfiguration *gridCfg, 
                 char absorber[],
                 double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                 double E_old_xdir[8][gridCfg->Ny][gridCfg->Nz], 
                 double E_old_ydir[gridCfg->Nx][8][gridCfg->Nz], 
                 double E_old_zdir[gridCfg->Nx][gridCfg->Ny][8] ) {
//{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk,
        offset;             // refers to EB_WAVE only

    double
        cnst;

    cnst    = (gridCfg->dt-gridCfg->dx)/(gridCfg->dt+gridCfg->dx);
    offset  = 2;

    // the string "absorber" is used to set which absorber is treated
    // the comparison is done with the strstr() function, which return the address
    // of the substring if found, NULL otherwise
    // NOTE: "if (strstr(absorber,"x1))" should be sufficient
    //       "if (strstr(absorber,"x1) != NULL)" should be equivalent

    // absorber into x-direction
    if ( strstr(absorber,"x1") ) {      
        //printf("abs_Mur_1st_v2: x1\n");
#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
                // absorber at x=0 grid boundary
                // Ex: odd-even-even
                EB_WAVE[offset+0+1][jj  ][kk  ] = E_old_xdir[2+1][jj  ][kk  ]
                    + cnst * (    EB_WAVE[offset+2+1][jj  ][kk  ] 
                              -E_old_xdir[0+1       ][jj  ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[offset+0  ][jj+1][kk  ] = E_old_xdir[2  ][jj+1][kk  ]
                    + cnst * (    EB_WAVE[offset+2  ][jj+1][kk  ] 
                              -E_old_xdir[0         ][jj+1][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[offset+0  ][jj  ][kk+1] = E_old_xdir[2  ][jj  ][kk+1]
                    + cnst * (    EB_WAVE[offset+2  ][jj  ][kk+1] 
                              -E_old_xdir[0         ][jj  ][kk+1] );
            }
        }
    }
    if ( strstr(absorber,"x2") ) {
        //printf("abs_Mur_1st_v2: x2\n");
#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
                // absorber at x=Nx grid boundary
                // Ex: odd-even-even
                EB_WAVE[gridCfg->Nx-2-offset+1][jj  ][kk  ]    = E_old_xdir[4+1][jj  ][kk  ]
                    + cnst * (    EB_WAVE[gridCfg->Nx-4-offset+1][jj  ][kk  ] 
                              -E_old_xdir[6+1                   ][jj  ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[gridCfg->Nx-2-offset  ][jj+1][kk  ]    = E_old_xdir[4  ][jj+1][kk  ]
                    + cnst * (    EB_WAVE[gridCfg->Nx-4-offset  ][jj+1][kk  ] 
                              -E_old_xdir[6                     ][jj+1][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[gridCfg->Nx-2-offset  ][jj  ][kk+1]    = E_old_xdir[4  ][jj  ][kk+1]
                    + cnst * (    EB_WAVE[gridCfg->Nx-4-offset  ][jj  ][kk+1] 
                              -E_old_xdir[6                     ][jj  ][kk+1] );
            }
        }
    }

    // absorber into y-direction
    if ( strstr(absorber,"y1") ) {
        //printf("abs_Mur_1st_v2: y1\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
        for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
            for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
                // absorber at y=0 grid boundary
                // Ex: odd-even-even
                EB_WAVE[ii+1][offset+0  ][kk  ] = E_old_ydir[ii+1][2  ][kk  ]
                    + cnst * (    EB_WAVE[ii+1][offset+2  ][kk  ]
                              -E_old_ydir[ii+1][0         ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[ii  ][offset+0+1][kk  ] = E_old_ydir[ii  ][2+1][kk  ]
                    + cnst * (    EB_WAVE[ii  ][offset+2+1][kk  ]
                              -E_old_ydir[ii  ][0+1       ][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[ii  ][offset+0  ][kk+1] = E_old_ydir[ii  ][2  ][kk+1]
                    + cnst * (    EB_WAVE[ii  ][offset+2  ][kk+1]
                              -E_old_ydir[ii  ][0         ][kk+1] );
            }
        }
    }
    if ( strstr(absorber,"y2") ) {
        //printf("abs_Mur_1st_v2: y2\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
        for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
            for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
                // absorber at y=Ny grid boundary
                // Ex: odd-even-even
                EB_WAVE[ii+1][gridCfg->Ny-2-offset  ][kk  ] = E_old_ydir[ii+1][4  ][kk  ]
                    + cnst * (    EB_WAVE[ii+1][gridCfg->Ny-4-offset  ][kk  ]
                              -E_old_ydir[ii+1][6                     ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[ii  ][gridCfg->Ny-2-offset+1][kk  ] = E_old_ydir[ii  ][4+1][kk  ]
                    + cnst * (    EB_WAVE[ii  ][gridCfg->Ny-4-offset+1][kk  ]
                              -E_old_ydir[ii  ][6+1                   ][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[ii  ][gridCfg->Ny-2-offset  ][kk+1] = E_old_ydir[ii  ][4  ][kk+1]
                    + cnst * (    EB_WAVE[ii  ][gridCfg->Ny-4-offset  ][kk+1]
                              -E_old_ydir[ii  ][6                     ][kk+1] );
            }
        }
    }

    // absorber into z-direction
    if ( strstr(absorber,"z1") ) {
        //printf("abs_Mur_1st_v2: z1\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
        for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
            for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
                // absorber at z=0 grid boundary
                // Ex: odd-even-even
                EB_WAVE[ii+1][jj  ][offset+0]   = E_old_zdir[ii+1][jj  ][2  ]
                    + cnst * (    EB_WAVE[ii+1][jj  ][offset+2  ]
                              -E_old_zdir[ii+1][jj  ][0  ]        );
                // Ey: even-odd-even
                EB_WAVE[ii  ][jj+1][offset+0]   = E_old_zdir[ii  ][jj+1][2  ]
                    + cnst * (    EB_WAVE[ii  ][jj+1][offset+2  ]
                              -E_old_zdir[ii  ][jj+1][0  ]        );
                // Ez: even-even-odd
                EB_WAVE[ii  ][jj  ][offset+0+1] = E_old_zdir[ii  ][jj  ][2+1]
                    + cnst * (    EB_WAVE[ii  ][jj  ][offset+2+1]
                              -E_old_zdir[ii  ][jj  ][0+1]        );
            }
        }
    }
    if ( strstr(absorber,"z2") ) {
        //printf("abs_Mur_1st_v2: z2\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
        for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
            for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
                // absorber at z=Nz grid boundary
                // Ex: odd-even-even
                EB_WAVE[ii+1][jj  ][gridCfg->Nz-2-offset  ]    = E_old_zdir[ii+1][jj  ][4  ]
                    + cnst * (    EB_WAVE[ii+1][jj  ][gridCfg->Nz-4-offset  ]
                              -E_old_zdir[ii+1][jj  ][6  ]     );
                // Ey: even-odd-even
                EB_WAVE[ii  ][jj+1][gridCfg->Nz-2-offset  ]    = E_old_zdir[ii  ][jj+1][4  ]
                    + cnst * (    EB_WAVE[ii  ][jj+1][gridCfg->Nz-4-offset  ]
                              -E_old_zdir[ii  ][jj+1][6  ]     );
                // Ez: even-even-odd
                EB_WAVE[ii  ][jj  ][gridCfg->Nz-2-offset+1]    = E_old_zdir[ii  ][jj  ][4+1]
                    + cnst * (    EB_WAVE[ii  ][jj  ][gridCfg->Nz-4-offset+1]
                              -E_old_zdir[ii  ][jj  ][6+1]     );
            }
        }
    }

    return EXIT_SUCCESS;

} //}}}


int abc_Mur_1st_ref( gridConfiguration *gridCfg,
                     double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                     double E_old_xdir[8][gridCfg->Ny][gridCfg->Nz_ref], 
                     double E_old_ydir[gridCfg->Nx][8][gridCfg->Nz_ref], 
                     double E_old_zdir[gridCfg->Nx][gridCfg->Ny][8] ) {
//{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk,
        offset;             // refers to EB_WAVE only

    double
        cnst;

    cnst    = (gridCfg->dt-gridCfg->dx)/(gridCfg->dt+gridCfg->dx);
    offset  = 2;

    // absorber into x-direction
#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            // absorber at x=0 grid boundary
            // Ex: odd-even-even
            EB_WAVE[offset+0+1][jj  ][kk  ] = E_old_xdir[2+1][jj  ][kk  ]
                + cnst * (    EB_WAVE[offset+2+1][jj  ][kk  ] 
                          -E_old_xdir[0+1       ][jj  ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[offset+0  ][jj+1][kk  ] = E_old_xdir[2  ][jj+1][kk  ]
                + cnst * (    EB_WAVE[offset+2  ][jj+1][kk  ] 
                          -E_old_xdir[0         ][jj+1][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[offset+0  ][jj  ][kk+1] = E_old_xdir[2  ][jj  ][kk+1]
                + cnst * (    EB_WAVE[offset+2  ][jj  ][kk+1] 
                          -E_old_xdir[0         ][jj  ][kk+1] );
            // absorber at x=Nx grid boundary
            // Ex: odd-even-even
            EB_WAVE[gridCfg->Nx-2-offset+1][jj  ][kk  ]    = E_old_xdir[4+1][jj  ][kk  ]
                + cnst * (    EB_WAVE[gridCfg->Nx-4-offset+1][jj  ][kk  ] 
                          -E_old_xdir[6+1                   ][jj  ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[gridCfg->Nx-2-offset  ][jj+1][kk  ]    = E_old_xdir[4  ][jj+1][kk  ]
                + cnst * (    EB_WAVE[gridCfg->Nx-4-offset  ][jj+1][kk  ] 
                          -E_old_xdir[6                     ][jj+1][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[gridCfg->Nx-2-offset  ][jj  ][kk+1]    = E_old_xdir[4  ][jj  ][kk+1]
                + cnst * (    EB_WAVE[gridCfg->Nx-4-offset  ][jj  ][kk+1] 
                          -E_old_xdir[6                     ][jj  ][kk+1] );
        }
    }

    // absorber into y-direction
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            // absorber at y=0 grid boundary
            // Ex: odd-even-even
            EB_WAVE[ii+1][offset+0  ][kk  ] = E_old_ydir[ii+1][2  ][kk  ]
                + cnst * (    EB_WAVE[ii+1][offset+2  ][kk  ]
                          -E_old_ydir[ii+1][0         ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[ii  ][offset+0+1][kk  ] = E_old_ydir[ii  ][2+1][kk  ]
                + cnst * (    EB_WAVE[ii  ][offset+2+1][kk  ]
                          -E_old_ydir[ii  ][0+1       ][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[ii  ][offset+0  ][kk+1] = E_old_ydir[ii  ][2  ][kk+1]
                + cnst * (    EB_WAVE[ii  ][offset+2  ][kk+1]
                          -E_old_ydir[ii  ][0         ][kk+1] );
            // absorber at y=Ny grid boundary
            // Ex: odd-even-even
            EB_WAVE[ii+1][gridCfg->Ny-2-offset  ][kk  ] = E_old_ydir[ii+1][4  ][kk  ]
                + cnst * (    EB_WAVE[ii+1][gridCfg->Ny-4-offset  ][kk  ]
                          -E_old_ydir[ii+1][6                     ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[ii  ][gridCfg->Ny-2-offset+1][kk  ] = E_old_ydir[ii  ][4+1][kk  ]
                + cnst * (    EB_WAVE[ii  ][gridCfg->Ny-4-offset+1][kk  ]
                          -E_old_ydir[ii  ][6+1                   ][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[ii  ][gridCfg->Ny-2-offset  ][kk+1] = E_old_ydir[ii  ][4  ][kk+1]
                + cnst * (    EB_WAVE[ii  ][gridCfg->Ny-4-offset  ][kk+1]
                          -E_old_ydir[ii  ][6                     ][kk+1] );
        }
    }

    // absorber into z-direction
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            // absorber at z=0 grid boundary
            // Ex: odd-even-even
            EB_WAVE[ii+1][jj  ][offset+0]   = E_old_zdir[ii+1][jj  ][2  ]
                + cnst * (    EB_WAVE[ii+1][jj  ][offset+2  ]
                          -E_old_zdir[ii+1][jj  ][0  ]        );
            // Ey: even-odd-even
            EB_WAVE[ii  ][jj+1][offset+0]   = E_old_zdir[ii  ][jj+1][2  ]
                + cnst * (    EB_WAVE[ii  ][jj+1][offset+2  ]
                          -E_old_zdir[ii  ][jj+1][0  ]        );
            // Ez: even-even-odd
            EB_WAVE[ii  ][jj  ][offset+0+1] = E_old_zdir[ii  ][jj  ][2+1]
                + cnst * (    EB_WAVE[ii  ][jj  ][offset+2+1]
                          -E_old_zdir[ii  ][jj  ][0+1]        );
            // absorber at z=Nz grid boundary
            // Ex: odd-even-even
            EB_WAVE[ii+1][jj  ][gridCfg->Nz_ref-2-offset  ]    = E_old_zdir[ii+1][jj  ][4  ]
                + cnst * (    EB_WAVE[ii+1][jj  ][gridCfg->Nz_ref-4-offset  ]
                          -E_old_zdir[ii+1][jj  ][6  ]     );
            // Ey: even-odd-even
            EB_WAVE[ii  ][jj+1][gridCfg->Nz_ref-2-offset  ]    = E_old_zdir[ii  ][jj+1][4  ]
                + cnst * (    EB_WAVE[ii  ][jj+1][gridCfg->Nz_ref-4-offset  ]
                          -E_old_zdir[ii  ][jj+1][6  ]     );
            // Ez: even-even-odd
            EB_WAVE[ii  ][jj  ][gridCfg->Nz_ref-2-offset+1]    = E_old_zdir[ii  ][jj  ][4+1]
                + cnst * (    EB_WAVE[ii  ][jj  ][gridCfg->Nz_ref-4-offset+1]
                          -E_old_zdir[ii  ][jj  ][6+1]     );
        }
    }

    return EXIT_SUCCESS;

} //}}}


double calc_poynt_4( gridConfiguration *gridCfg, 
                     int pwr_dect, char absorber[],
                     double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                     double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref] ) {
//{{{

    size_t
        ii, jj, kk;
    double
        poynt;

    poynt   = .0;

    // P = E x H
    // Px = Ey*Hz - Ez*Hy
    // Py = Ez*Hx - Ex*Hz
    // Pz = Ex*Hy - Ey*Hx
    // Bx: even-odd-odd
    // By: odd-even-odd
    // Bz: odd-odd-even
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd
    
    if ( strcmp(absorber,"ref_z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(gridCfg->Nx-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<=(gridCfg->Ny-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE_ref[ii+1][jj  ][pwr_dect  ]
                          *EB_WAVE_ref[ii+1][jj  ][pwr_dect+1]
                          -EB_WAVE_ref[ii  ][jj+1][pwr_dect  ]
                          *EB_WAVE_ref[ii  ][jj+1][pwr_dect+1] );
            }
        }
    } else if ( strcmp(absorber,"z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(gridCfg->Nx-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<=(gridCfg->Ny-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( ( EB_WAVE[ii+1][jj  ][pwr_dect  ] - EB_WAVE_ref[ii+1][jj  ][pwr_dect  ] )
                          *( EB_WAVE[ii+1][jj  ][pwr_dect+1] - EB_WAVE_ref[ii+1][jj  ][pwr_dect+1] )
                          -( EB_WAVE[ii  ][jj+1][pwr_dect  ] - EB_WAVE_ref[ii  ][jj+1][pwr_dect  ] ) 
                          *( EB_WAVE[ii  ][jj+1][pwr_dect+1] - EB_WAVE_ref[ii  ][jj+1][pwr_dect+1] ) );
            }
        }
    } else if ( strcmp(absorber,"z2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(gridCfg->Nx-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<=(gridCfg->Ny-pwr_dect-2) ; jj+=2) {
                // z2-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE[ii+1][jj  ][gridCfg->Nz-pwr_dect  ]
                          *EB_WAVE[ii+1][jj  ][gridCfg->Nz-pwr_dect-1]
                          -EB_WAVE[ii  ][jj+1][gridCfg->Nz-pwr_dect  ]
                          *EB_WAVE[ii  ][jj+1][gridCfg->Nz-pwr_dect-1] );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(gridCfg->Ny-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(gridCfg->Nz-pwr_dect-2) ; kk+=2) {
                // x1-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[pwr_dect  ][jj+1][kk  ]
                          *EB_WAVE[pwr_dect+1][jj+1][kk  ]
                          -EB_WAVE[pwr_dect  ][jj  ][kk+1]
                          *EB_WAVE[pwr_dect+1][jj  ][kk+1] );
            }
        }
    } else if ( strcmp(absorber,"x2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(gridCfg->Ny-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(gridCfg->Nz-pwr_dect-2) ; kk+=2) {
                // x2-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[gridCfg->Nx-pwr_dect  ][jj+1][kk  ]
                          *EB_WAVE[gridCfg->Nx-pwr_dect-1][jj+1][kk  ]
                          -EB_WAVE[gridCfg->Nx-pwr_dect  ][jj  ][kk+1]
                          *EB_WAVE[gridCfg->Nx-pwr_dect-1][jj  ][kk+1] );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(gridCfg->Nx-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(gridCfg->Nz-pwr_dect-2) ; kk+=2) {
                // y1-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][pwr_dect  ][kk+1]
                          *EB_WAVE[ii  ][pwr_dect+1][kk+1]
                          -EB_WAVE[ii+1][pwr_dect  ][kk  ]
                          *EB_WAVE[ii+1][pwr_dect+1][kk  ] );
            }
        }
    } else if ( strcmp(absorber,"y2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(gridCfg->Nx-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(gridCfg->Nz-pwr_dect-2) ; kk+=2) {
                // y2-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][gridCfg->Ny-pwr_dect  ][kk+1]
                          *EB_WAVE[ii  ][gridCfg->Ny-pwr_dect-1][kk+1]
                          -EB_WAVE[ii+1][gridCfg->Ny-pwr_dect  ][kk  ]
                          *EB_WAVE[ii+1][gridCfg->Ny-pwr_dect-1][kk  ] );
            }
        }
    }
    
    return fabs(poynt);
} //}}}


double calc_poynt_5( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] ) {
//{{{
// same as calc_poynt_4, but with spatial averaging of B-component

    size_t
        ii, jj, kk;
    double
        poynt;

    poynt   = .0;

    // P = E x H
    // Px = Ey*Hz - Ez*Hy
    // Py = Ez*Hx - Ex*Hz
    // Pz = Ex*Hy - Ey*Hx
    // Bx: even-odd-odd
    // By: odd-even-odd
    // Bz: odd-odd-even
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd
    
    if ( strcmp(absorber,"ref_z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE_ref[ii+1][jj  ][pwr_dect  ]
                          *( EB_WAVE_ref[ii-1][jj  ][pwr_dect+1]
                            +EB_WAVE_ref[ii+3][jj  ][pwr_dect+1]
                            +EB_WAVE_ref[ii+1][jj-2][pwr_dect+1]
                            +EB_WAVE_ref[ii+1][jj+2][pwr_dect+1] )*.25
                          -EB_WAVE_ref[ii  ][jj+1][pwr_dect  ]
                          *( EB_WAVE_ref[ii-2][jj+1][pwr_dect+1]
                            +EB_WAVE_ref[ii+2][jj+1][pwr_dect+1]
                            +EB_WAVE_ref[ii  ][jj-1][pwr_dect+1]
                            +EB_WAVE_ref[ii  ][jj+3][pwr_dect+1] )*.25 );
            }
        }
    } else if ( strcmp(absorber,"z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( ( EB_WAVE[ii+1][jj  ][pwr_dect  ] - EB_WAVE_ref[ii+1][jj  ][pwr_dect  ] )
                          *( ( EB_WAVE[ii-1][jj  ][pwr_dect+1] 
                              +EB_WAVE[ii+3][jj  ][pwr_dect+1] 
                              +EB_WAVE[ii+1][jj-2][pwr_dect+1] 
                              +EB_WAVE[ii+1][jj+2][pwr_dect+1] )*.25
                            -( EB_WAVE_ref[ii-1][jj  ][pwr_dect+1]
                              +EB_WAVE_ref[ii+3][jj  ][pwr_dect+1]
                              +EB_WAVE_ref[ii+1][jj-2][pwr_dect+1]
                              +EB_WAVE_ref[ii+1][jj+2][pwr_dect+1] )*.25 )
                          -( EB_WAVE[ii  ][jj+1][pwr_dect  ] - EB_WAVE_ref[ii  ][jj+1][pwr_dect  ] ) 
                          *( ( EB_WAVE[ii-2][jj+1][pwr_dect+1] 
                              +EB_WAVE[ii+2][jj+1][pwr_dect+1] 
                              +EB_WAVE[ii  ][jj-1][pwr_dect+1] 
                              +EB_WAVE[ii  ][jj+3][pwr_dect+1] )*.25
                            -( EB_WAVE_ref[ii-2][jj+1][pwr_dect+1]
                              +EB_WAVE_ref[ii+2][jj+1][pwr_dect+1]
                              +EB_WAVE_ref[ii  ][jj-1][pwr_dect+1]
                              +EB_WAVE_ref[ii  ][jj+3][pwr_dect+1] )*.25 ) );
            }
        }
    } else if ( strcmp(absorber,"z2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z2-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE[ii+1][jj  ][N_z-pwr_dect  ]
                          *( EB_WAVE[ii-1][jj  ][N_z-pwr_dect-1]
                            +EB_WAVE[ii+3][jj  ][N_z-pwr_dect-1]
                            +EB_WAVE[ii+1][jj-2][N_z-pwr_dect-1]
                            +EB_WAVE[ii+1][jj+2][N_z-pwr_dect-1] )*.25
                          -EB_WAVE[ii  ][jj+1][N_z-pwr_dect  ]
                          *( EB_WAVE[ii-2][jj+1][N_z-pwr_dect-1]
                            +EB_WAVE[ii+2][jj+1][N_z-pwr_dect-1]
                            +EB_WAVE[ii  ][jj-1][N_z-pwr_dect-1]
                            +EB_WAVE[ii  ][jj+3][N_z-pwr_dect-1] )*.25 );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x1-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[pwr_dect  ][jj+1][kk  ]
                          *( EB_WAVE[pwr_dect+1][jj-1][kk  ]
                            +EB_WAVE[pwr_dect+1][jj+3][kk  ]
                            +EB_WAVE[pwr_dect+1][jj+1][kk-2]
                            +EB_WAVE[pwr_dect+1][jj+1][kk+2] )*.25
                          -EB_WAVE[pwr_dect  ][jj  ][kk+1]
                          *( EB_WAVE[pwr_dect+1][jj-2][kk+1]
                            +EB_WAVE[pwr_dect+1][jj+2][kk+1]
                            +EB_WAVE[pwr_dect+1][jj  ][kk-1]
                            +EB_WAVE[pwr_dect+1][jj  ][kk+3] )*.25 );
            }
        }
    } else if ( strcmp(absorber,"x2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x2-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[N_x-pwr_dect  ][jj+1][kk  ]
                          *( EB_WAVE[N_x-pwr_dect-1][jj-1][kk  ]
                            +EB_WAVE[N_x-pwr_dect-1][jj+3][kk  ]
                            +EB_WAVE[N_x-pwr_dect-1][jj+1][kk-2]
                            +EB_WAVE[N_x-pwr_dect-1][jj+1][kk+2] )*.25
                          -EB_WAVE[N_x-pwr_dect  ][jj  ][kk+1]
                          *( EB_WAVE[N_x-pwr_dect-1][jj-2][kk+1]
                            +EB_WAVE[N_x-pwr_dect-1][jj+2][kk+1]
                            +EB_WAVE[N_x-pwr_dect-1][jj  ][kk-1]
                            +EB_WAVE[N_x-pwr_dect-1][jj  ][kk+3] )*.25 );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y1-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][pwr_dect  ][kk+1]
                          *( EB_WAVE[ii-2][pwr_dect+1][kk+1]
                            +EB_WAVE[ii+2][pwr_dect+1][kk+1]
                            +EB_WAVE[ii  ][pwr_dect+1][kk-1]
                            +EB_WAVE[ii  ][pwr_dect+1][kk+3] )*.25
                          -EB_WAVE[ii+1][pwr_dect  ][kk  ]
                          *( EB_WAVE[ii-1][pwr_dect+1][kk  ]
                            +EB_WAVE[ii+3][pwr_dect+1][kk  ]
                            +EB_WAVE[ii+1][pwr_dect+1][kk-2]
                            +EB_WAVE[ii+1][pwr_dect+1][kk+2] )*.25 );
            }
        }
    } else if ( strcmp(absorber,"y2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y2-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][N_y-pwr_dect  ][kk+1]
                          *( EB_WAVE[ii-2][N_y-pwr_dect-1][kk+1]
                            +EB_WAVE[ii+2][N_y-pwr_dect-1][kk+1]
                            +EB_WAVE[ii  ][N_y-pwr_dect-1][kk-1]
                            +EB_WAVE[ii  ][N_y-pwr_dect-1][kk+3] )*.25
                          -EB_WAVE[ii+1][N_y-pwr_dect  ][kk  ]
                          *( EB_WAVE[ii-1][N_y-pwr_dect-1][kk  ]
                            +EB_WAVE[ii+3][N_y-pwr_dect-1][kk  ]
                            +EB_WAVE[ii+1][N_y-pwr_dect-1][kk-2]
                            +EB_WAVE[ii+1][N_y-pwr_dect-1][kk+2] )*.25 );
            }
        }
    }
    
    return fabs(poynt);
} //}}}


double calc_poynt_6( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] ) {
//{{{
// same as calc_poynt_5, but full spatial averaging of B-component

    size_t
        ii, jj, kk;
    double
        poynt;

    poynt   = .0;

    // P = E x H
    // Px = Ey*Hz - Ez*Hy
    // Py = Ez*Hx - Ex*Hz
    // Pz = Ex*Hy - Ey*Hx
    // Bx: even-odd-odd
    // By: odd-even-odd
    // Bz: odd-odd-even
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd
    
    if ( strcmp(absorber,"ref_z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE_ref[ii+1][jj  ][pwr_dect  ]
                          *( EB_WAVE_ref[ii-1][jj  ][pwr_dect+1]
                            +EB_WAVE_ref[ii+3][jj  ][pwr_dect+1]
                            +EB_WAVE_ref[ii+1][jj-2][pwr_dect+1]
                            +EB_WAVE_ref[ii+1][jj+2][pwr_dect+1] 
                            +EB_WAVE_ref[ii+1][jj  ][pwr_dect-1]
                            +EB_WAVE_ref[ii+1][jj  ][pwr_dect+3] )/6.
                          -EB_WAVE_ref[ii  ][jj+1][pwr_dect  ]
                          *( EB_WAVE_ref[ii-2][jj+1][pwr_dect+1]
                            +EB_WAVE_ref[ii+2][jj+1][pwr_dect+1]
                            +EB_WAVE_ref[ii  ][jj-1][pwr_dect+1]
                            +EB_WAVE_ref[ii  ][jj+3][pwr_dect+1]
                            +EB_WAVE_ref[ii  ][jj+1][pwr_dect-1]
                            +EB_WAVE_ref[ii  ][jj+1][pwr_dect+3] )/6. );
            }
        }
    } else if ( strcmp(absorber,"z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( ( EB_WAVE[ii+1][jj  ][pwr_dect  ] - EB_WAVE_ref[ii+1][jj  ][pwr_dect  ] )
                          *( ( EB_WAVE[ii-1][jj  ][pwr_dect+1] 
                              +EB_WAVE[ii+3][jj  ][pwr_dect+1] 
                              +EB_WAVE[ii+1][jj-2][pwr_dect+1] 
                              +EB_WAVE[ii+1][jj+2][pwr_dect+1]
                              +EB_WAVE[ii+1][jj  ][pwr_dect-1] 
                              +EB_WAVE[ii+1][jj  ][pwr_dect+3] )/6.
                            -( EB_WAVE_ref[ii-1][jj  ][pwr_dect+1]
                              +EB_WAVE_ref[ii+3][jj  ][pwr_dect+1]
                              +EB_WAVE_ref[ii+1][jj-2][pwr_dect+1]
                              +EB_WAVE_ref[ii+1][jj+2][pwr_dect+1]
                              +EB_WAVE_ref[ii+1][jj  ][pwr_dect-1]
                              +EB_WAVE_ref[ii+1][jj  ][pwr_dect+3] )/6. )
                          -( EB_WAVE[ii  ][jj+1][pwr_dect  ] - EB_WAVE_ref[ii  ][jj+1][pwr_dect  ] ) 
                          *( ( EB_WAVE[ii-2][jj+1][pwr_dect+1] 
                              +EB_WAVE[ii+2][jj+1][pwr_dect+1] 
                              +EB_WAVE[ii  ][jj-1][pwr_dect+1] 
                              +EB_WAVE[ii  ][jj+3][pwr_dect+1]
                              +EB_WAVE[ii  ][jj+1][pwr_dect-1] 
                              +EB_WAVE[ii  ][jj+1][pwr_dect+3] )/6.
                            -( EB_WAVE_ref[ii-2][jj+1][pwr_dect+1]
                              +EB_WAVE_ref[ii+2][jj+1][pwr_dect+1]
                              +EB_WAVE_ref[ii  ][jj-1][pwr_dect+1]
                              +EB_WAVE_ref[ii  ][jj+3][pwr_dect+1]
                              +EB_WAVE_ref[ii  ][jj+1][pwr_dect-1]
                              +EB_WAVE_ref[ii  ][jj+1][pwr_dect+3] )/6. ) );
            }
        }
    } else if ( strcmp(absorber,"z2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z2-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE[ii+1][jj  ][N_z-pwr_dect  ]
                          *( EB_WAVE[ii-1][jj  ][N_z-pwr_dect-1]
                            +EB_WAVE[ii+3][jj  ][N_z-pwr_dect-1]
                            +EB_WAVE[ii+1][jj-2][N_z-pwr_dect-1]
                            +EB_WAVE[ii+1][jj+2][N_z-pwr_dect-1]
                            +EB_WAVE[ii+1][jj  ][N_z-pwr_dect-3]
                            +EB_WAVE[ii+1][jj  ][N_z-pwr_dect+1] )/6.
                          -EB_WAVE[ii  ][jj+1][N_z-pwr_dect  ]
                          *( EB_WAVE[ii-2][jj+1][N_z-pwr_dect-1]
                            +EB_WAVE[ii+2][jj+1][N_z-pwr_dect-1]
                            +EB_WAVE[ii  ][jj-1][N_z-pwr_dect-1]
                            +EB_WAVE[ii  ][jj+3][N_z-pwr_dect-1]
                            +EB_WAVE[ii  ][jj+1][N_z-pwr_dect-3]
                            +EB_WAVE[ii  ][jj+1][N_z-pwr_dect+1] )/6. );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x1-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[pwr_dect  ][jj+1][kk  ]
                          *( EB_WAVE[pwr_dect+1][jj-1][kk  ]
                            +EB_WAVE[pwr_dect+1][jj+3][kk  ]
                            +EB_WAVE[pwr_dect+1][jj+1][kk-2]
                            +EB_WAVE[pwr_dect+1][jj+1][kk+2]
                            +EB_WAVE[pwr_dect-1][jj+1][kk  ]
                            +EB_WAVE[pwr_dect+3][jj+1][kk  ] )/6.
                          -EB_WAVE[pwr_dect  ][jj  ][kk+1]
                          *( EB_WAVE[pwr_dect+1][jj-2][kk+1]
                            +EB_WAVE[pwr_dect+1][jj+2][kk+1]
                            +EB_WAVE[pwr_dect+1][jj  ][kk-1]
                            +EB_WAVE[pwr_dect+1][jj  ][kk+3]
                            +EB_WAVE[pwr_dect-1][jj  ][kk+1]
                            +EB_WAVE[pwr_dect+3][jj  ][kk+1] )/6. );
            }
        }
    } else if ( strcmp(absorber,"x2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x2-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[N_x-pwr_dect  ][jj+1][kk  ]
                          *( EB_WAVE[N_x-pwr_dect-1][jj-1][kk  ]
                            +EB_WAVE[N_x-pwr_dect-1][jj+3][kk  ]
                            +EB_WAVE[N_x-pwr_dect-1][jj+1][kk-2]
                            +EB_WAVE[N_x-pwr_dect-1][jj+1][kk+2]
                            +EB_WAVE[N_x-pwr_dect-3][jj+1][kk  ]
                            +EB_WAVE[N_x-pwr_dect+1][jj+1][kk  ] )/6.
                          -EB_WAVE[N_x-pwr_dect  ][jj  ][kk+1]
                          *( EB_WAVE[N_x-pwr_dect-1][jj-2][kk+1]
                            +EB_WAVE[N_x-pwr_dect-1][jj+2][kk+1]
                            +EB_WAVE[N_x-pwr_dect-1][jj  ][kk-1]
                            +EB_WAVE[N_x-pwr_dect-1][jj  ][kk+3]
                            +EB_WAVE[N_x-pwr_dect-3][jj  ][kk+1]
                            +EB_WAVE[N_x-pwr_dect+1][jj  ][kk+1] )/6. );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y1-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][pwr_dect  ][kk+1]
                          *( EB_WAVE[ii-2][pwr_dect+1][kk+1]
                            +EB_WAVE[ii+2][pwr_dect+1][kk+1]
                            +EB_WAVE[ii  ][pwr_dect+1][kk-1]
                            +EB_WAVE[ii  ][pwr_dect+1][kk+3]
                            +EB_WAVE[ii  ][pwr_dect-1][kk+1]
                            +EB_WAVE[ii  ][pwr_dect+3][kk+1] )/6.
                          -EB_WAVE[ii+1][pwr_dect  ][kk  ]
                          *( EB_WAVE[ii-1][pwr_dect+1][kk  ]
                            +EB_WAVE[ii+3][pwr_dect+1][kk  ]
                            +EB_WAVE[ii+1][pwr_dect+1][kk-2]
                            +EB_WAVE[ii+1][pwr_dect+1][kk+2]
                            +EB_WAVE[ii+1][pwr_dect-1][kk  ]
                            +EB_WAVE[ii+1][pwr_dect+3][kk  ] )/6. );
            }
        }
    } else if ( strcmp(absorber,"y2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y2-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][N_y-pwr_dect  ][kk+1]
                          *( EB_WAVE[ii-2][N_y-pwr_dect-1][kk+1]
                            +EB_WAVE[ii+2][N_y-pwr_dect-1][kk+1]
                            +EB_WAVE[ii  ][N_y-pwr_dect-1][kk-1]
                            +EB_WAVE[ii  ][N_y-pwr_dect-1][kk+3]
                            +EB_WAVE[ii  ][N_y-pwr_dect-3][kk+1]
                            +EB_WAVE[ii  ][N_y-pwr_dect+1][kk+1] )/6.
                          -EB_WAVE[ii+1][N_y-pwr_dect  ][kk  ]
                          *( EB_WAVE[ii-1][N_y-pwr_dect-1][kk  ]
                            +EB_WAVE[ii+3][N_y-pwr_dect-1][kk  ]
                            +EB_WAVE[ii+1][N_y-pwr_dect-1][kk-2]
                            +EB_WAVE[ii+1][N_y-pwr_dect-1][kk+2]
                            +EB_WAVE[ii+1][N_y-pwr_dect-3][kk  ]
                            +EB_WAVE[ii+1][N_y-pwr_dect+1][kk  ] )/6. );
            }
        }
    }
    
    return fabs(poynt);
} //}}}


double calc_poynt_7( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] ) {
//{{{
// same as calc_poynt_4, but with symmetrizing of B-component

    size_t
        ii, jj, kk;
    double
        poynt;

    poynt   = .0;

    // P = E x H
    // Px = Ey*Hz - Ez*Hy
    // Py = Ez*Hx - Ex*Hz
    // Pz = Ex*Hy - Ey*Hx
    // Bx: even-odd-odd
    // By: odd-even-odd
    // Bz: odd-odd-even
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd
    
    if ( strcmp(absorber,"ref_z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE_ref[ii+1][jj  ][pwr_dect  ]
                          *( EB_WAVE_ref[ii+1][jj  ][pwr_dect+1]
                            +EB_WAVE_ref[ii+1][jj  ][pwr_dect-1] )*.5
                          -EB_WAVE_ref[ii  ][jj+1][pwr_dect  ]
                          *( EB_WAVE_ref[ii  ][jj+1][pwr_dect+1]
                            +EB_WAVE_ref[ii  ][jj+1][pwr_dect-1] )*.5 );
            }
        }
    } else if ( strcmp(absorber,"z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( ( EB_WAVE[ii+1][jj  ][pwr_dect  ] - EB_WAVE_ref[ii+1][jj  ][pwr_dect  ] )
                          *( ( EB_WAVE[ii+1][jj  ][pwr_dect+1] 
                              +EB_WAVE[ii+1][jj  ][pwr_dect-1] )*.5
                            -( EB_WAVE_ref[ii+1][jj  ][pwr_dect+1]
                              +EB_WAVE_ref[ii+1][jj  ][pwr_dect-1] )*.5 )
                          -( EB_WAVE[ii  ][jj+1][pwr_dect  ] - EB_WAVE_ref[ii  ][jj+1][pwr_dect  ] ) 
                          *( ( EB_WAVE[ii  ][jj+1][pwr_dect+1] 
                              +EB_WAVE[ii  ][jj+1][pwr_dect-1] )*.5
                            -( EB_WAVE_ref[ii  ][jj+1][pwr_dect+1]
                              +EB_WAVE_ref[ii  ][jj+1][pwr_dect-1] )*.5 ) );
            }
        }
    } else if ( strcmp(absorber,"z2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z2-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE[ii+1][jj  ][N_z-pwr_dect  ]
                          *( EB_WAVE[ii+1][jj  ][N_z-pwr_dect-1]
                            +EB_WAVE[ii+1][jj  ][N_z-pwr_dect+1] )*.5
                          -EB_WAVE[ii  ][jj+1][N_z-pwr_dect  ]
                          *( EB_WAVE[ii  ][jj+1][N_z-pwr_dect-1]
                            +EB_WAVE[ii  ][jj+1][N_z-pwr_dect+1] )*.5 );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x1-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[pwr_dect  ][jj+1][kk  ]
                          *( EB_WAVE[pwr_dect+1][jj+1][kk  ]
                            +EB_WAVE[pwr_dect-1][jj+1][kk  ] )*.5
                          -EB_WAVE[pwr_dect  ][jj  ][kk+1]
                          *( EB_WAVE[pwr_dect+1][jj  ][kk+1] 
                            +EB_WAVE[pwr_dect-1][jj  ][kk+1] )*.5 );
            }
        }
    } else if ( strcmp(absorber,"x2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x2-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[N_x-pwr_dect  ][jj+1][kk  ]
                          *( EB_WAVE[N_x-pwr_dect-1][jj+1][kk  ]
                            +EB_WAVE[N_x-pwr_dect+1][jj+1][kk  ] )*.5
                          -EB_WAVE[N_x-pwr_dect  ][jj  ][kk+1]
                          *( EB_WAVE[N_x-pwr_dect-1][jj  ][kk+1] 
                            +EB_WAVE[N_x-pwr_dect+1][jj  ][kk+1] )*.5 );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y1-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][pwr_dect  ][kk+1]
                          *( EB_WAVE[ii  ][pwr_dect+1][kk+1]
                            +EB_WAVE[ii  ][pwr_dect-1][kk+1] )*.5
                          -EB_WAVE[ii+1][pwr_dect  ][kk  ]
                          *( EB_WAVE[ii+1][pwr_dect+1][kk  ]
                            +EB_WAVE[ii+1][pwr_dect-1][kk  ] )*.5 );
            }
        }
    } else if ( strcmp(absorber,"y2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y2-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][N_y-pwr_dect  ][kk+1]
                          *( EB_WAVE[ii  ][N_y-pwr_dect-1][kk+1]
                            +EB_WAVE[ii  ][N_y-pwr_dect+1][kk+1] )*.5
                          -EB_WAVE[ii+1][N_y-pwr_dect  ][kk  ]
                          *( EB_WAVE[ii+1][N_y-pwr_dect-1][kk  ]
                            +EB_WAVE[ii+1][N_y-pwr_dect+1][kk  ] )*.5 );
            }
        }
    }
    
    return fabs(poynt);
} //}}}


double calc_power_EE_1( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                       int pwr_dect, char absorber[],
                       double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] ) {
//{{{
    // sum-up the power (squared E-field) in the absorbers

    size_t
        ii, jj, kk;
    double
        power;

    power   = .0;

    // Bx: even-odd-odd
    // By: odd-even-odd
    // Bz: odd-odd-even
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    if ( strcmp(absorber,"ref_z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:power)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                power += ( pow(EB_WAVE_ref[ii+1][jj  ][pwr_dect  ], 2)
                          +pow(EB_WAVE_ref[ii  ][jj+1][pwr_dect  ], 2)
                          +pow(EB_WAVE_ref[ii  ][jj  ][pwr_dect+1], 2) );
            }
        }
    } else if ( strcmp(absorber,"z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:power)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                power += ( pow((EB_WAVE[ii+1][jj  ][pwr_dect  ] - EB_WAVE_ref[ii+1][jj  ][pwr_dect  ]), 2)
                          +pow((EB_WAVE[ii  ][jj+1][pwr_dect  ] - EB_WAVE_ref[ii  ][jj+1][pwr_dect  ]), 2)
                          +pow((EB_WAVE[ii  ][jj  ][pwr_dect+1] - EB_WAVE_ref[ii  ][jj  ][pwr_dect+1]), 2) );
            }
        }
    } else if ( strcmp(absorber,"z2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:power)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z2-plane
                power += ( pow(EB_WAVE[ii+1][jj  ][N_z-pwr_dect  ], 2)
                          +pow(EB_WAVE[ii  ][jj+1][N_z-pwr_dect  ], 2)
                          +pow(EB_WAVE[ii  ][jj  ][N_z-pwr_dect-1], 2) );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:power)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x1-plane
                power += ( pow(EB_WAVE[pwr_dect+1][jj  ][kk  ], 2)
                          +pow(EB_WAVE[pwr_dect  ][jj+1][kk  ], 2)
                          +pow(EB_WAVE[pwr_dect  ][jj  ][kk+1], 2) );
            }
        }
    } else if ( strcmp(absorber,"x2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:power)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x2-plane
                power += ( pow(EB_WAVE[N_x-pwr_dect-1][jj  ][kk  ], 2)
                          +pow(EB_WAVE[N_x-pwr_dect  ][jj+1][kk  ], 2)
                          +pow(EB_WAVE[N_x-pwr_dect  ][jj  ][kk+1], 2) );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:power)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y1-plane
                power += ( pow(EB_WAVE[ii+1][pwr_dect  ][kk  ], 2)
                          +pow(EB_WAVE[ii  ][pwr_dect+1][kk  ], 2)
                          +pow(EB_WAVE[ii  ][pwr_dect  ][kk+1], 2) );
            }
        }
    } else if ( strcmp(absorber,"y2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:power)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y2-plane
                power += ( pow(EB_WAVE[ii+1][N_y-pwr_dect  ][kk  ], 2)
                          +pow(EB_WAVE[ii  ][N_y-pwr_dect+1][kk  ], 2)
                          +pow(EB_WAVE[ii  ][N_y-pwr_dect  ][kk+1], 2) );
            }
        }
    }
    
    return fabs(power);
} //}}}


#ifdef DETECTOR_ANTENNA_1D
int detAnt1D_storeValues( gridConfiguration *gridCfg, 
                          size_t detAnt_ypos, size_t detAnt_zpos,
                          int tt, 
                          double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                          double detAnt_fields[gridCfg->Nx/2][5] ) { 
    //{{{
    size_t
        ii;

    double
        foo;

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

#pragma omp parallel default(shared) private(ii,foo)
#pragma omp for
    for ( ii=2 ; ii <= gridCfg->Nx-2 ; ii+=2 ) {
        // calculate abs(E)
        foo = sqrt(  pow(EB_WAVE[ii+1][detAnt_ypos  ][detAnt_zpos  ],2)
                    +pow(EB_WAVE[ii  ][detAnt_ypos+1][detAnt_zpos  ],2)
                    +pow(EB_WAVE[ii  ][detAnt_ypos  ][detAnt_zpos+1],2) );

        // sum of E over time
        // Ex*Ex
        detAnt_fields[ii/2][0]  += pow( EB_WAVE[ii+1][detAnt_ypos  ][detAnt_zpos  ], 2 );
        // Ey*Ey
        detAnt_fields[ii/2][1]  += pow( EB_WAVE[ii  ][detAnt_ypos+1][detAnt_zpos  ], 2 );
        // Ez*Ez
        detAnt_fields[ii/2][2]  += pow( EB_WAVE[ii  ][detAnt_ypos  ][detAnt_zpos+1], 2 );
        // E*E
        detAnt_fields[ii/2][3]  += foo*foo;

        // corresponding to an rms(E)-like quantity
        detAnt_fields[ii/2][4]  += ( foo * sqrt(1./( (double)(tt)/(double)(gridCfg->period) + 1e-6 )) );

        //printf( "tt = %d, ii = %d, sum_t(E*E) = %13.5e\n",
        //        tt, ii, detAnt_fields[ii/2][3] );
    }

    return EXIT_SUCCESS;

}//}}}
#endif


#if defined(HDF5) && defined(DETECTOR_ANTENNA_1D)
int detAnt1D_write2hdf5( int N_x, 
                         char filename[], char detAnt_groupName[], 
                         size_t detAnt_ypos, size_t detAnt_zpos,
                         double detAnt_fields[N_x/2][5] ){
    //#{{{

    // hdf related variables
    hid_t       file_id, dataset_id, dataspace_id,      // object identifiers
                group_id__detAnt,
                dataspace_id_i, 
                dataset_id_i;
    hsize_t     dims[1];                                // size used for dimensions
    herr_t      status;                                 // function return value

    // hdf5 related variables for applying shuffle and gzip filter
    hid_t           dcpl;
    hsize_t         chunk[1];
    unsigned int    filter_info;
    int             filter_avail;

    // required for check if hdf5-file already exists
    struct      stat st;

    int         ii;

    double       data2save[N_x/2];

    dims[0]     = N_x/2;
    chunk[0]    = N_x/2;

    // assume as default setting that filters are available
    filter_avail = 1;

    //data2save = dvector( 0, n_elem/2 );
    set2zero_1D( N_x/2, data2save );

    printf( "Will write data for detector antenna position y=%05ld, z=%05ld into file %s\n", 
            detAnt_ypos, detAnt_zpos, filename );

    // check if specified hdf5 file already exists
    // if not, create new one; if yes, open and add dataset to it
    if ( stat( filename, &st )==0 ) {
        // open file for read + write access
        file_id = H5Fopen( filename,            // filename
                           H5F_ACC_RDWR,        // allow read & write access (_RDONLY for read only)
                           H5P_DEFAULT);        // file access property list (default one)
    } else {
        // create a new file using default properties.
        file_id = H5Fcreate( filename,          // filename
                             H5F_ACC_TRUNC,     // how file should be created (removes existing file)
                             H5P_DEFAULT,       // file creating property list
                             H5P_DEFAULT);      // file access property list
    }

    // create group for different data to be stored
    group_id__detAnt = H5Gcreate2( file_id,
                                   detAnt_groupName,
                                   H5P_DEFAULT,
                                   H5P_DEFAULT,
                                   H5P_DEFAULT);

    // check if filters for shuffle and gzip exists, if yes, apply them
    // (check is done, because the filters are optional parts of hdf5 library)
    // check if gzip-filter is available
    if ( !(H5Zfilter_avail( H5Z_FILTER_DEFLATE )) ) {
        printf( "WARNING: gzip filter not available (for hdf5)\n" );
        filter_avail = 0;
    } else {
        status = H5Zget_filter_info( H5Z_FILTER_DEFLATE, &filter_info );
        if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ||
             !(filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED) ) {
            printf( "WARNING: gzip filter not available for encoding and decoding (for hdf5)\n" );
            filter_avail = 0;
        }
        if (status < 0) printf( "ERROR: could not get hdf5 filter info\n" );
    }
    // check if shuffle-filter is available
    if ( !(H5Zfilter_avail( H5Z_FILTER_SHUFFLE )) ) {
        printf( "WARNING: shuffle filter not available (for hdf5)\n" );
        filter_avail = 0;
    } else {
        status = H5Zget_filter_info( H5Z_FILTER_SHUFFLE, &filter_info );
        if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ||
             !(filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED) ) {
            printf( "WARNING: shuffle filter not available for encoding and decoding (for hdf5)\n" );
            filter_avail = 0;
        }
        if (status < 0) printf( "ERROR: could not get hdf5 filter info\n" );
    }

    // apply shuffle and gzip filters, if available
    if (filter_avail) {
        // create dataset creation property list
        dcpl = H5Pcreate( H5P_DATASET_CREATE );
        // add shuffle filter and gzip compression filter
        // note that the order of filter is significant: first shuffle!
        // order of filters applied correspond to order in which they are invoked when writin gdata
        status = H5Pset_shuffle( dcpl );
        if (status < 0) printf( "ERROR: could not apply shuffle filter\n" );
        status = H5Pset_deflate( dcpl, 9 );
        if (status < 0) printf( "ERROR: could not apply gzip filter\n" );
        status = H5Pset_chunk(dcpl, 1, chunk );
        if (status < 0) printf( "ERROR: could not set size of chunk\n" );
    } 

    // store spatial coordinate
    // prepare array
    for ( ii=2 ; ii<=N_x-2 ; ii+=2 ) {
        data2save[ii/2] = (double)(ii) ;
    }
    // create data space
    dataspace_id_i = H5Screate_simple( 1,       // rank of array (number of dimensions of dataspace)
                                       dims,    // array of the size of each dimension
                                       NULL);   // allow stretching of data space (NULL=no)
//    printf( "dataspace_id_i=%d\n", dataspace_id_i);
    // create new dataset and links it to location in file
    printf( "start to create dataset 'j'\n" );
    if (filter_avail) {
        dataset_id_i = H5Dcreate( group_id__detAnt,     // file identifier (or group identifier)
                                  "i",                  // name of dataset (relative to group specified, if speficied)
                                  H5T_NATIVE_DOUBLE,    // datatype to use when creating dataset
                                  dataspace_id_i,       // dataspace identifier
                                  H5P_DEFAULT,          // link creation property list
                                  dcpl,                 // dataset creation property list
                                  H5P_DEFAULT);         // dataset access property list
    } else {
        dataset_id_i = H5Dcreate( group_id__detAnt,     // file identifier (or group identifier)
                                  "i",                  // name of dataset (relative to group specified, if speficied)
                                  H5T_NATIVE_DOUBLE,    // datatype to use when creating dataset
                                  dataspace_id_i,       // dataspace identifier
                                  H5P_DEFAULT,
                                  H5P_DEFAULT,
                                  H5P_DEFAULT);
    }
    // write the dataset
    status = H5Dwrite( dataset_id_i,        // dataset identifier
                       H5T_NATIVE_DOUBLE,   // informs hdf about format of data in memory of computer
                       H5S_ALL,             // identifier of memory dataspace
                       H5S_ALL,             // file space identifier
                       H5P_DEFAULT,         // data transfer property list
                       data2save);          // pointer to data array
    if (status < 0) printf( "ERROR: could not write dataset 'i'\n" );

    status       = H5Dclose(dataset_id_i);
    if (status < 0) printf( "ERROR: could not close dataset 'i'\n" );
    status       = H5Sclose(dataspace_id_i);
    if (status < 0) printf( "ERROR: could not close dataspace for 'i' dataset\n" );


    // store position
    dims[0] = 1;
    data2save[0] = (double)(detAnt_ypos);
    data2save[1] = (double)(detAnt_zpos);
    dataspace_id = H5Screate_simple( 1, dims, NULL); 
    // detAnt_ypos
    printf( "start to create dataset 'detAnt_ypos'\n" );
    dataset_id   = H5Dcreate( group_id__detAnt, "detAnt_ypos", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
    status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data2save[0]);
    if (status < 0) printf( "ERROR: could not write dataset 'detAnt_ypos'\n" );
    status       = H5Dclose(dataset_id);
    if (status < 0) printf( "ERROR: could not close dataset 'detAnt_ypos'\n" );
    // detAnt_zpos
    printf( "start to create dataset 'detAnt_zpos'\n" );
    dataset_id   = H5Dcreate( group_id__detAnt, "detAnt_zpos", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
    status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data2save[1]);
    if (status < 0) printf( "ERROR: could not write dataset 'detAnt_zpos'\n" );
    status       = H5Dclose(dataset_id);
    if (status < 0) printf( "ERROR: could not close dataset 'detAnt_zpos'\n" );
    status       = H5Sclose(dataspace_id);
    if (status < 0) printf( "ERROR: could not close dataspace for detAnt-position datasets\n" );

    // store sum_ExEx
    dims[0] = N_x/2;
    // since all following arrays have same dimension, dataspace_id needs to be created only once
    // and not closed with H5Sclose(dataspace_id) after each dataset
    dataspace_id = H5Screate_simple( 1, dims, NULL);
    // prepare array to be saved
    set2zero_1D( N_x/2, data2save );
    for ( ii=2 ; ii<=N_x-2 ; ii+=2 )
        data2save[ii/2] = detAnt_fields[ii/2][0];
    printf( "start to create dataset 'sum_ExEx'\n" );
    if (filter_avail)
        dataset_id   = H5Dcreate( group_id__detAnt, "sum_ExEx", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, dcpl, H5P_DEFAULT);  
    else
        dataset_id   = H5Dcreate( group_id__detAnt, "sum_ExEx", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
    status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2save);
    if (status < 0) printf( "ERROR: could not write dataset 'sum_ExEx'\n" );
    status       = H5Dclose(dataset_id);
    if (status < 0) printf( "ERROR: could not close dataset 'sum_ExEx'\n" );

    // store sum_EyEy 
    set2zero_1D( N_x/2, data2save );
    for ( ii=2 ; ii<=N_x-2 ; ii+=2 )
        data2save[ii/2] = detAnt_fields[ii/2][1];
    printf( "start to create dataset 'sum_EyEy'\n" );
    if (filter_avail)
        dataset_id   = H5Dcreate( group_id__detAnt, "sum_EyEy", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, dcpl, H5P_DEFAULT);  
    else
        dataset_id   = H5Dcreate( group_id__detAnt, "sum_EyEy", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
    status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2save);
    if (status < 0) printf( "ERROR: could not write dataset 'sum_EyEy'\n" );
    status       = H5Dclose(dataset_id);
    if (status < 0) printf( "ERROR: could not close dataset 'sum_EyEy'\n" );

    // store sum_EzEz
    set2zero_1D( N_x/2, data2save );
    for ( ii=2 ; ii<=N_x-2 ; ii+=2 )
        data2save[ii/2] = detAnt_fields[ii/2][2];
    printf( "start to create dataset 'sum_EzEz'\n" );
    if (filter_avail)
        dataset_id   = H5Dcreate( group_id__detAnt, "sum_EzEz", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, dcpl, H5P_DEFAULT);  
    else
        dataset_id   = H5Dcreate( group_id__detAnt, "sum_EzEz", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
    status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2save);
    if (status < 0) printf( "ERROR: could not write dataset 'sum_EzEz'\n" );
    status       = H5Dclose(dataset_id);
    if (status < 0) printf( "ERROR: could not close dataset 'sum_EzEz'\n" );
    
    // store sum_EE
    set2zero_1D( N_x/2, data2save );
    for ( ii=2 ; ii<=N_x-2 ; ii+=2 )
        data2save[ii/2] = detAnt_fields[ii/2][3];
    printf( "start to create dataset 'sum_EE'\n" );
    if (filter_avail)
        dataset_id   = H5Dcreate( group_id__detAnt, "sum_EE", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, dcpl, H5P_DEFAULT);  
    else
        dataset_id   = H5Dcreate( group_id__detAnt, "sum_EE", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
    status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2save);
    if (status < 0) printf( "ERROR: could not write dataset 'sum_EE'\n" );
    status       = H5Dclose(dataset_id);
    if (status < 0) printf( "ERROR: could not close dataset 'sum_EE'\n" );

    // store rmsE
    set2zero_1D( N_x/2, data2save );
    for ( ii=2 ; ii<=N_x-2 ; ii+=2 )
        data2save[ii/2] = detAnt_fields[ii/2][4];
    printf( "start to create dataset 'rms_E'\n" );
    if (filter_avail)
        dataset_id   = H5Dcreate( group_id__detAnt, "rms_E", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, dcpl, H5P_DEFAULT);  
    else
        dataset_id   = H5Dcreate( group_id__detAnt, "rms_E", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
    status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2save);
    if (status < 0) printf( "ERROR: could not write dataset 'rmsE'\n" );
    status       = H5Dclose(dataset_id);
    if (status < 0) printf( "ERROR: could not close dataset 'rmsE'\n" );

    status       = H5Sclose(dataspace_id);
    if (status < 0) printf( "ERROR: could not close dataspace for datasets of E-fields\n" );
    
    // terminate access and free ressources/identifiers
    if (filter_avail) {
        status = H5Pclose( dcpl );
        if (status < 0) printf( "ERROR: could not close filter\n" );
    }
    status = H5Gclose( group_id__detAnt );
    if (status < 0) printf( "ERROR: could not close group detAnt\n" );
    // file 
    status = H5Fclose(file_id);
    if (status < 0) printf( "ERROR: could not close group file '%s'\n", filename );

    return EXIT_SUCCESS;
}//#}}}
#endif


int set2zero_1D( size_t N_x, double arr_1D[N_x] ){
//{{{

    size_t
        ii;

#pragma omp parallel for default(shared) private(ii)
    for (ii=0 ; ii<N_x ; ++ii) {
        arr_1D[ii] = .0;
    }

    return EXIT_SUCCESS;
} //}}}


int set2zero_3D( size_t N_x, size_t N_y, size_t N_z, double arr_3D[N_x][N_y][N_z] ){
//{{{

    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii<N_x ; ++ii) {
        for (jj=0 ; jj<N_y ; ++jj) {
            for (kk=0 ; kk<N_z ; ++kk) {
                arr_3D[ii][jj][kk]  = .0;
            }
        }
    }

    return EXIT_SUCCESS;
} //}}}


