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

// define structures
typedef struct gridConfiguration {
    int
        Nx, Ny, Nz, 
        Nz_ref,
        d_absorb,
        t_end;
    double
        period,
        dx,dt;
} gridConfiguration;
typedef struct beamConfiguration {
    int
        ant_x, ant_y, ant_z;
    double
        antAngle_zy, antAngle_zx,
        ant_w0x, ant_w0y;
} beamConfiguration;

// prototyping
int make_antenna_profile( gridConfiguration *gridCfg, beamConfiguration *beamCfg,
                          int exc_signal, 
                          double z2waist,
                          double antField_xy[gridCfg->Nx/2][gridCfg->Ny/2], double antPhaseTerms[gridCfg->Nx/2][gridCfg->Ny/2] );
int make_density_profile( gridConfiguration *gridCfg, 
                          int ne_profile, 
                          double cntrl_para, 
                          double n_e[gridCfg->Nx/2][gridCfg->Ny/2][gridCfg->Nz/2] );
int set_densityInAbsorber_v2( double period, int d_absorb, 
                              char absorber[], 
                              size_t N_x, size_t N_y, size_t N_z,
                              double n_e[N_x/2][N_y/2][N_z/2] );
int make_B0_profile( gridConfiguration *gridCfg,
                     int B0_profile, 
                     double cntrl_para, 
                     double J_B0[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] );
int add_source( beamConfiguration *beamCfg, 
                int exc_signal,
                size_t N_x, size_t N_y, size_t N_z,  
                double period, 
                int t_int, double omega_t, 
                double antField_xy[N_x/2][N_y/2], double antPhaseTerms[N_x/2][N_y/2],
                double EB_WAVE[N_x][N_y][N_z] );
int apply_absorber( gridConfiguration *gridCfg, 
                    double eco, 
                    double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] );
int apply_absorber_ref( gridConfiguration *gridCfg, 
                        double eco, 
                        double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref] );
int apply_absorber_v2( size_t N_x, size_t N_y, size_t N_z, int d_absorb, double eco, 
                       char absorber[],
                       double EB_WAVE[N_x][N_y][N_z] );
int abc_Mur_saveOldE_xdir( size_t N_x, size_t N_y, size_t N_z, 
                           double EB_WAVE[N_x][N_y][N_z], 
                           double E_old[8][N_y][N_z] );
int abc_Mur_saveOldE_ydir( size_t N_x, size_t N_y, size_t N_z, 
                           double EB_WAVE[N_x][N_y][N_z], 
                           double E_old[N_y][8][N_z] );
int abc_Mur_saveOldE_zdir( size_t N_x, size_t N_y, size_t N_z, 
                           double EB_WAVE[N_x][N_y][N_z], 
                           double E_old[N_x][N_y][8] );
int abs_Mur_1st( size_t N_x, size_t N_y, size_t N_z,
                 double dt, double dx, 
                 double EB_WAVE[N_x][N_y][N_z], 
                 double E_old_xdir[8][N_y][N_z], double E_old_ydir[N_x][8][N_z], double E_old_zdir[N_x][N_y][8] );
int abs_Mur_1st_v2( size_t N_x, size_t N_y, size_t N_z,
                    double dt, double dx, 
                    char absorber[],
                    double EB_WAVE[N_x][N_y][N_z], 
                    double E_old_xdir[8][N_y][N_z], double E_old_ydir[N_x][8][N_z], double E_old_zdir[N_x][N_y][8] );
int advance_J( size_t dim1, size_t dim2, size_t dim3, 
               double EB_WAVE[dim1][dim2][dim3], double J_B0[dim1][dim2][dim3],
               double n_e[dim1/2][dim2/2][dim3/2], 
               double dt );
int advance_B( size_t dim1, size_t dim2, size_t dim3, 
               double EB_WAVE[dim1][dim2][dim3], 
               double dx, double dt );
int advance_E( size_t dim1, size_t dim2, size_t dim3, 
               double EB_WAVE[dim1][dim2][dim3], double J_B0[dim1][dim2][dim3],
               double dx, double dt );
int advance_E_vacuum( size_t dim1, size_t dim2, size_t dim3, 
                      double EB_WAVE[dim1][dim2][dim3],
                      double dx, double dt );
int calc_poynt_1( size_t N_x, size_t N_y, size_t N_z, 
                  int pwr_dect, char absorber[], double poynt[3],  
                  double EB_WAVE[N_x][N_y][N_z] );
double calc_poynt_2( size_t N_x, size_t N_y, size_t N_z, 
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z] );
double calc_poynt_3( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] );
double calc_poynt_4( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] );
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
int writeTimetraces2ascii( int dim0, int dim1, int t_end, double period, 
                           char filename[], double timetraces[dim0][dim1] );
#ifdef DETECTOR_ANTENNA_1D
int detAnt1D_storeValues( size_t N_x, size_t N_y, size_t N_z,
                          size_t detAnt_ypos, size_t detAnt_zpos,
                          int tt, double period,  
                          double EB_WAVE[N_x][N_y][N_z], double detAnt_fields[N_x/2][5] );
#endif
#if defined(HDF5) && defined(DETECTOR_ANTENNA_1D)
int detAnt1D_write2hdf5( int N_x, 
                         char filename[], char detAnt_groupName[], 
                         size_t detAnt_ypos, size_t detAnt_zpos,
                         double detAnt_fields[N_x/2][5] );
#endif
#ifdef HDF5
int writeMyHDF_v4( int dim0, int dim1, int dim2, char filename[], char dataset[], double array_3D[dim0][dim1][dim2] );
int writeConfig2HDF( char filename[], int N_x, int N_y, int N_z, int period, int d_absorb );
int readMyHDF( int dim0, int dim1, int dim2, char filename[], char dataset[], double array_3D[dim0][dim1][dim2]);
#endif


int main( int argc, char *argv[] ) {
//{{{

    struct gridConfiguration gridCfg;
    struct beamConfiguration beamCfg;

    int
        ii,jj,kk,
        t_int, t_end, T_wave, 

        scale,
        NX, NY, NZ, 
        NZ_ref,

#ifdef _OPENMP
        n_threads,                          // number of threads that will be used (OpenMP)
#endif
        d_absorb, pwr_dect,

#ifdef DETECTOR_ANTENNA_1D
        detAnt_01_zpos,
        detAnt_02_zpos,
        detAnt_03_zpos,
        detAnt_04_zpos,
        detAnt_05_zpos,
        detAnt_06_zpos,
        detAnt_07_zpos,
        detAnt_01_ypos,
#endif

        len_str,                            // return value of snprintf
        opt_ret;                            // return value of getopt (reading input parameter)

    double
        period,
        dx,dt,
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

        aux,

        omega_t;

    char
        dSet_name[PATH_MAX],
        filename_hdf5[PATH_MAX];            // filename of hdf5 file for output

    bool
        angle_zx_set,                       // is antAngle_zx set during call ?
        angle_zy_set;                       // is antAngle_zy set during call ?

    // set-up grid
    scale           = 1;
    period          = 16*scale;
    gridCfg.period  = 16*scale;
#if BOUNDARY == 1
    d_absorb        = (int)(3*period);
    gridCfg.d_absorb= (int)(3*period);
#elif BOUNDARY == 2
    d_absorb        = 8;
    gridCfg.d_absorb= 8;
#endif
    NX          = (280)*scale;
    NY          = (220)*scale;
    NZ          = (160)*scale;
    gridCfg.Nx  = (280)*scale;
    gridCfg.Ny  = (220)*scale;
    gridCfg.Nz  = (160)*scale;
    NZ_ref          = 2*d_absorb + (int)period;
    gridCfg.Nz_ref  = 2*d_absorb + (int)period;
    t_end           = (int)((30)*period);
    gridCfg.t_end   = (int)((30)*period);

    // arrays realized as variable-length array (VLA)
    // E- and B-wavefield
    double (*EB_WAVE)[NY][NZ]           = calloc(NX, sizeof *EB_WAVE);
    double (*EB_WAVE_ref)[NY][NZ_ref]   = calloc(NX, sizeof *EB_WAVE_ref);
    // J-wavefield (in plasma) and background magnetic field
    double (*J_B0)[NY][NZ]              = calloc(NX, sizeof *J_B0);
    // background electron plasma density
    double (*n_e)[NY/2][NZ/2]           = calloc(NX/2, sizeof *n_e);
    // used when writing data into hdf5-files
    double (*data2save)[NY/2][NZ/2]     = calloc(NX/2, sizeof *data2save);
    // antenna: envelope of injected field
    double (*antField_xy)[NY/2]         = calloc(NX/2, sizeof *antField_xy);
    // antenna: phase terms 
    double (*antPhaseTerms)[NY/2]       = calloc(NX/2, sizeof *antPhaseTerms);
    // time traces
    double (*timetraces)[8]             = calloc((t_end/(int)period), sizeof *timetraces);

    // old E-fields required for Mur's boundary condition
#if BOUNDARY == 2
    double (*E_Xdir_OLD)[NY][NZ]        = calloc(8,  sizeof *E_Xdir_OLD);
    double (*E_Ydir_OLD)[8][NZ]         = calloc(NX, sizeof *E_Ydir_OLD);
    double (*E_Zdir_OLD)[NY][8]         = calloc(NX, sizeof *E_Zdir_OLD);
    double (*E_Xdir_OLD_ref)[NY][NZ_ref]= calloc(8,  sizeof *E_Xdir_OLD_ref);
    double (*E_Ydir_OLD_ref)[8][NZ_ref] = calloc(NX, sizeof *E_Ydir_OLD_ref);
    double (*E_Zdir_OLD_ref)[NY][8]     = calloc(NX, sizeof *E_Zdir_OLD_ref);
#endif

    // array for detector antennas
    // sum_t(Ex*Ex) | sum_t(Ey*Ey) | sum_t(Ez*Ez) | sum_t(E*E) | rms(E)
#ifdef DETECTOR_ANTENNA_1D
    // TODO: change into 3D array, such that each detector antenna corresponds
    //       to one 2D array; that way it can be written much more failsafe...
    //       requires some changes in procedures for storing and saving
    double (*detAnt_01_fields)[5]       = calloc(NX, sizeof *detAnt_01_fields);
    double (*detAnt_02_fields)[5]       = calloc(NX, sizeof *detAnt_02_fields);
    double (*detAnt_03_fields)[5]       = calloc(NX, sizeof *detAnt_03_fields);
    double (*detAnt_04_fields)[5]       = calloc(NX, sizeof *detAnt_04_fields);
    double (*detAnt_05_fields)[5]       = calloc(NX, sizeof *detAnt_05_fields);
    double (*detAnt_06_fields)[5]       = calloc(NX, sizeof *detAnt_06_fields);
    double (*detAnt_07_fields)[5]       = calloc(NX, sizeof *detAnt_07_fields);
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


    beamCfg.ant_x       = NX/2;
    beamCfg.ant_y       = NY/2;
    beamCfg.ant_z       = d_absorb + 4;
    // positions have to be even numbers, to ensure fields are accessed correctly
    if ((beamCfg.ant_x % 2) != 0)  ++beamCfg.ant_x;
    if ((beamCfg.ant_y % 2) != 0)  ++beamCfg.ant_y;
    if ((beamCfg.ant_z % 2) != 0)  ++beamCfg.ant_z;
    beamCfg.ant_w0x     = 2;
    beamCfg.ant_w0y     = 2;

    pwr_dect    = d_absorb;

#ifdef DETECTOR_ANTENNA_1D
    detAnt_01_ypos  = beamCfg.ant_y;
    detAnt_01_zpos  = beamCfg.ant_z+2;
    detAnt_02_zpos  = round(beamCfg.ant_z+2 + 1*4.67*period); // steps of 5 cm for 28 GHz
    detAnt_03_zpos  = round(beamCfg.ant_z+2 + 2*4.67*period);
    detAnt_04_zpos  = round(beamCfg.ant_z+2 + 3*4.67*period);
    detAnt_05_zpos  = round(beamCfg.ant_z+2 + 4*4.67*period);
    detAnt_06_zpos  = round(beamCfg.ant_z+2 + 5*4.67*period);
    detAnt_07_zpos  = round(beamCfg.ant_z+2 + 6*4.67*period);
    // positions have to be even numbers, to ensure fields are accessed correctly
    if ((detAnt_01_ypos % 2) != 0)  ++detAnt_01_ypos;
    if ((detAnt_01_zpos % 2) != 0)  ++detAnt_01_zpos;
    if ((detAnt_02_zpos % 2) != 0)  ++detAnt_02_zpos;
    if ((detAnt_03_zpos % 2) != 0)  ++detAnt_03_zpos;
    if ((detAnt_04_zpos % 2) != 0)  ++detAnt_04_zpos;
    if ((detAnt_05_zpos % 2) != 0)  ++detAnt_05_zpos;
    if ((detAnt_06_zpos % 2) != 0)  ++detAnt_06_zpos;
    if ((detAnt_07_zpos % 2) != 0)  ++detAnt_07_zpos;
    // issue a warning when detector antenna position is beyond Nz
    if (detAnt_07_zpos > (NZ - d_absorb)) {
        printf( "ERROR: check the detector antenna positions into z direction\n" );
        printf( "       NZ-d_absorb = %d, detAnt_07_zpos = %d", 
                NZ-d_absorb, detAnt_07_zpos );
    }
#endif

    // dt/dx = 0.5 is commenly used in 2D FDTD codes
    // Note that period refers to the wavelength in the numerical grid and not
    // in the "physical" grid (where one grid cell is equal to one Yee cell).
    // This means that in the physical grid, the wavelength is period/2, thus
    // in the equations we have to use period/2 for the wavelength.
    dx      = 1./(period/2);
    dt      = 1./(2.*(period/2));
        
#if BOUNDARY == 1
    eco         = 10./(double)(period);
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
                            1, 
                          -(298.87)*.0,                // .2/l_0*period = -298.87
                          antField_xy, antPhaseTerms );
    printf( "...done defining antenna field\n" );

    printf( "starting defining background plasma density\n" );
    make_density_profile( &gridCfg,  
            // ne_profile: 1 = plasma mirror
            //             2 = linearly increasing profile
            2,
            // cntrl_para: ne_profile=1 --> 0: plane mirror; oblique mirror: -.36397; 20 degrees: -.17633
            //             ne_profile=2 --> k0*Ln: 25
            25,
            n_e );
    printf( " ...setting density in absorber to 0...\n ");
    //set_densityInAbsorber_v2( period, d_absorb, "z1", NX, NY, NZ, n_e );
    //set_densityInAbsorber_v2( period, d_absorb, "x1x2y1y2z1", NX, NY, NZ, n_e );
    printf( "...done defining background plasma density\n" );

    printf( "starting defining background magnetic field...\n" );
    // B0x: even-odd-odd
    // B0y: odd-even-odd
    // B0z: odd-odd-even
    make_B0_profile(
            &gridCfg,
            // B0_profile: 1 = constant field
            1, 
            // cntrl_para: B0_profile=1 --> value of Y
            .85, 
            J_B0 );
    printf( "...done defining background magnetic field\n" );

    // print some info to console
    printf( "Nx = %d, Ny = %d, Nz = %d\n", NX, NY, NZ );
    printf( "period = %d\n", (int)(period) );
    printf( "d_absorb = %d\n", d_absorb );
    printf( "t_end = %d\n", (int)(t_end) );
    printf( "antAngle_zx = %.2f, antAngle_zy = %.2f\n", beamCfg.antAngle_zx, beamCfg.antAngle_zy );
    printf( "ant_w0x = %.2f, ant_w0y = %.2f\n", beamCfg.ant_w0x, beamCfg.ant_w0y ); 
    printf( "ant_x = %d, ant_y = %d, ant_z = %d\n", beamCfg.ant_x, beamCfg.ant_y, beamCfg.ant_z );
    printf( "Boundary condition set to '%d'\n", BOUNDARY );
#ifdef DETECTOR_ANTENNA_1D
    printf( "detector antenna positions: z1 = %d, y1 = %d\n", detAnt_01_zpos, detAnt_01_ypos );
    printf( "detector antenna positions: z2 = %d, y1 = %d\n", detAnt_02_zpos, detAnt_01_ypos );
    printf( "detector antenna positions: z3 = %d, y1 = %d\n", detAnt_03_zpos, detAnt_01_ypos );
    printf( "detector antenna positions: z4 = %d, y1 = %d\n", detAnt_04_zpos, detAnt_01_ypos );
    printf( "detector antenna positions: z5 = %d, y1 = %d\n", detAnt_05_zpos, detAnt_01_ypos );
    printf( "detector antenna positions: z6 = %d, y1 = %d\n", detAnt_06_zpos, detAnt_01_ypos );
    printf( "detector antenna positions: z7 = %d, y1 = %d\n", detAnt_07_zpos, detAnt_01_ypos );
#endif

#ifdef _OPENMP
#pragma omp parallel private(n_threads)
    {
    n_threads = omp_get_num_threads();
    printf( "number of threads that will be used (OpenMP) = %d\n", n_threads );
    }
#endif


    for (t_int=0 ; t_int <=t_end ; ++t_int) {
        
        omega_t += 2.*M_PI/period;

        // to avoid precision problems when a lot of pi's are summed up        
        if (omega_t >= 2.*M_PI) {
            omega_t    += -2.*M_PI;
            T_wave     += 1;
            //printf("status: number of oscillation periods: %d (t_int= %d) \n",T_wave,t_int);
        }

        // add source
        add_source( &beamCfg,
                    3,
        //add_source( 1,
                    NX, NY, NZ, period, t_int, omega_t, 
                    antField_xy, antPhaseTerms, EB_WAVE );
        add_source( &beamCfg,
                    3,
        //add_source( 1,
                    NX, NY, NZ_ref, period, t_int, omega_t, 
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
        advance_J( NX, NY, NZ, EB_WAVE, J_B0, n_e, dt );

        // advance B
        advance_B( NX, NY, NZ,     EB_WAVE,     dx, dt );
        advance_B( NX, NY, NZ_ref, EB_WAVE_ref, dx, dt );
        
        // advance E
        advance_E(        NX, NY, NZ,     EB_WAVE,     J_B0, dx, dt );
        advance_E_vacuum( NX, NY, NZ_ref, EB_WAVE_ref,       dx, dt );

        // apply Mur's boundary conditions
#if BOUNDARY == 2
        abs_Mur_1st_v2( NX, NY, NZ, dt, dx, "x1x2y1y2z1z2",  
                        EB_WAVE, E_Xdir_OLD, E_Ydir_OLD, E_Zdir_OLD );
        abs_Mur_1st( NX, NY, NZ_ref, dt, dx, 
                     EB_WAVE_ref, E_Xdir_OLD_ref, E_Ydir_OLD_ref, E_Zdir_OLD_ref );
        abc_Mur_saveOldE_xdir( NX, NY, NZ, EB_WAVE, E_Xdir_OLD );
        abc_Mur_saveOldE_ydir( NX, NY, NZ, EB_WAVE, E_Ydir_OLD );
        abc_Mur_saveOldE_zdir( NX, NY, NZ, EB_WAVE, E_Zdir_OLD );
        abc_Mur_saveOldE_xdir( NX, NY, NZ_ref, EB_WAVE_ref, E_Xdir_OLD_ref );
        abc_Mur_saveOldE_ydir( NX, NY, NZ_ref, EB_WAVE_ref, E_Ydir_OLD_ref );
        abc_Mur_saveOldE_zdir( NX, NY, NZ_ref, EB_WAVE_ref, E_Zdir_OLD_ref );
#endif

#ifdef DETECTOR_ANTENNA_1D
        // store wavefields for detector antennas over the final 10 
        // oscillation periods, it was found previously that only one period
        // does not result in a too nice average
        if ( t_int >= (t_end-10*period) ) {
            detAnt1D_storeValues( NX, NY, NZ, detAnt_01_ypos, detAnt_01_zpos,
                                  t_int, period,  
                                  EB_WAVE, detAnt_01_fields );
            detAnt1D_storeValues( NX, NY, NZ, detAnt_01_ypos, detAnt_02_zpos,
                                  t_int, period,  
                                  EB_WAVE, detAnt_02_fields );
            detAnt1D_storeValues( NX, NY, NZ, detAnt_01_ypos, detAnt_03_zpos,
                                  t_int, period,  
                                  EB_WAVE, detAnt_03_fields );
            detAnt1D_storeValues( NX, NY, NZ, detAnt_01_ypos, detAnt_04_zpos,
                                  t_int, period,  
                                  EB_WAVE, detAnt_04_fields );
            detAnt1D_storeValues( NX, NY, NZ, detAnt_01_ypos, detAnt_05_zpos,
                                  t_int, period,  
                                  EB_WAVE, detAnt_05_fields );
            detAnt1D_storeValues( NX, NY, NZ, detAnt_01_ypos, detAnt_06_zpos,
                                  t_int, period,  
                                  EB_WAVE, detAnt_06_fields );
            detAnt1D_storeValues( NX, NY, NZ, detAnt_01_ypos, detAnt_07_zpos,
                                  t_int, period,  
                                  EB_WAVE, detAnt_07_fields );
        }
#endif

        // IQ detector for power detection
        if ( t_int >= 20*period ) {
            // z1-plane and z2-plane
            poynt_z1_ref    = calc_poynt_4( NX, NY, NZ, NZ_ref, pwr_dect, "ref_z1", EB_WAVE, EB_WAVE_ref );
            poynt_z1        = calc_poynt_4( NX, NY, NZ, NZ_ref, pwr_dect, "z1",     EB_WAVE, EB_WAVE_ref );
            //poynt_z1_ref    = ( calc_poynt_4( NX, NY, NZ, NZ_ref, pwr_dect  , "ref_z1", EB_WAVE, EB_WAVE_ref )
            //                   +calc_poynt_4( NX, NY, NZ, NZ_ref, pwr_dect+2, "ref_z1", EB_WAVE, EB_WAVE_ref ) )*.5;
            //poynt_z1        = ( calc_poynt_4( NX, NY, NZ, NZ_ref, pwr_dect,   "z1",     EB_WAVE, EB_WAVE_ref )
            //                   +calc_poynt_4( NX, NY, NZ, NZ_ref, pwr_dect+2, "z1",     EB_WAVE, EB_WAVE_ref ) )*.5;
            poynt_z2        = calc_poynt_4( NX, NY, NZ, NZ_ref, pwr_dect, "z2",     EB_WAVE, EB_WAVE_ref );
            // x1-plane and x2-plane
            //poynt_x1        = calc_poynt_4( NX, NY, NZ, NZ_ref, pwr_dect, "x1", EB_WAVE, EB_WAVE_ref );
            //poynt_x2        = calc_poynt_4( NX, NY, NZ, NZ_ref, pwr_dect, "x2", EB_WAVE, EB_WAVE_ref );
            // y1-plane and y2-plane
            //poynt_y1        = calc_poynt_4( NX, NY, NZ, NZ_ref, pwr_dect, "y1", EB_WAVE, EB_WAVE_ref );
            //poynt_y2        = calc_poynt_4( NX, NY, NZ, NZ_ref, pwr_dect, "y2", EB_WAVE, EB_WAVE_ref );

            
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
            power_EE_ref    += calc_power_EE_1( NX, NY, NZ, NZ_ref, d_absorb, "ref_z1", EB_WAVE, EB_WAVE_ref );
            power_EE_z1     += calc_power_EE_1( NX, NY, NZ, NZ_ref, d_absorb, "z1",     EB_WAVE, EB_WAVE_ref );
            power_EE_z2     += calc_power_EE_1( NX, NY, NZ, NZ_ref, d_absorb, "z2",     EB_WAVE, EB_WAVE_ref );
            // x1-plane and x2-plane
            power_EE_x1     += calc_power_EE_1( NX, NY, NZ, NZ_ref, d_absorb, "x1",     EB_WAVE, EB_WAVE_ref );
            power_EE_x2     += calc_power_EE_1( NX, NY, NZ, NZ_ref, d_absorb, "x2",     EB_WAVE, EB_WAVE_ref );
            // y1-plane and y2-plane
            power_EE_y1     += calc_power_EE_1( NX, NY, NZ, NZ_ref, d_absorb, "y1",     EB_WAVE, EB_WAVE_ref );
            power_EE_y2     += calc_power_EE_1( NX, NY, NZ, NZ_ref, d_absorb, "y2",     EB_WAVE, EB_WAVE_ref );
            */

        }


        if ( (t_int % (int)(period)) == 4 )  {
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
    for ( ii=0 ; ii<(t_end/(int)period) ; ++ii )
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
    writeTimetraces2ascii( (t_end/(int)period), 8, t_end, period, 
                           "timetraces2.dat", timetraces );

    // save into hdf5
    // abs(E)
    // prepare array for that
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii<NX ; ii+=2) {
        for (jj=0 ; jj<NY ; jj+=2) {
            for (kk=0 ; kk<NZ ; kk+=2) {
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
        printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( NX/2, NY/2, NZ/2, filename_hdf5, dSet_name, data2save) ) ;
    }
    // density
    sprintf( dSet_name, "n_e" );
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( NX/2, NY/2, NZ/2, filename_hdf5, dSet_name, n_e) ) ;
    set2zero_3D( NX/2, NY/2, NZ/2, data2save );
    // background magnetic field
    // B0x: even-odd-odd
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii<NX ; ii+=2) {
        for (jj=0 ; jj<NY ; jj+=2) {
            for (kk=0 ; kk<NZ ; kk+=2) {
                data2save[(ii/2)][(jj/2)][(kk/2)] = J_B0[ii  ][jj+1][kk+1];
            }
        }
    }
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( NX/2, NY/2, NZ/2, filename_hdf5, "B0x", data2save) ) ;
    set2zero_3D( NX/2, NY/2, NZ/2, data2save );
    // B0y: odd-even-odd
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii<NX ; ii+=2) {
        for (jj=0 ; jj<NY ; jj+=2) {
            for (kk=0 ; kk<NZ ; kk+=2) {
                data2save[(ii/2)][(jj/2)][(kk/2)] = J_B0[ii+1][jj  ][kk+1];
            }
        }
    }
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( NX/2, NY/2, NZ/2, filename_hdf5, "B0y", data2save) ) ;
    set2zero_3D( NX/2, NY/2, NZ/2, data2save );
    // B0z: odd-odd-even
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii<NX ; ii+=2) {
        for (jj=0 ; jj<NY ; jj+=2) {
            for (kk=0 ; kk<NZ ; kk+=2) {
                data2save[(ii/2)][(jj/2)][(kk/2)] = J_B0[ii+1][jj+1][kk  ];
            }
        }
    }
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( NX/2, NY/2, NZ/2, filename_hdf5, "B0z", data2save) ) ;
    set2zero_3D( NX/2, NY/2, NZ/2, data2save );

    writeConfig2HDF( filename_hdf5, NX, NY, NZ, period, d_absorb );


#if defined(HDF5) && defined(DETECTOR_ANTENNA_1D)
    detAnt1D_write2hdf5( NX, filename_hdf5, "/detAnt_01" , 
                         detAnt_01_ypos, detAnt_01_zpos,
                         detAnt_01_fields );
    detAnt1D_write2hdf5( NX, filename_hdf5, "/detAnt_02" , 
                         detAnt_01_ypos, detAnt_02_zpos,
                         detAnt_02_fields );
    detAnt1D_write2hdf5( NX, filename_hdf5, "/detAnt_03" , 
                         detAnt_01_ypos, detAnt_03_zpos,
                         detAnt_03_fields );
    detAnt1D_write2hdf5( NX, filename_hdf5, "/detAnt_04" , 
                         detAnt_01_ypos, detAnt_04_zpos,
                         detAnt_04_fields );
    detAnt1D_write2hdf5( NX, filename_hdf5, "/detAnt_05" , 
                         detAnt_01_ypos, detAnt_05_zpos,
                         detAnt_05_fields );
    detAnt1D_write2hdf5( NX, filename_hdf5, "/detAnt_06" , 
                         detAnt_01_ypos, detAnt_06_zpos,
                         detAnt_06_fields );
    detAnt1D_write2hdf5( NX, filename_hdf5, "/detAnt_07" , 
                         detAnt_01_ypos, detAnt_07_zpos,
                         detAnt_07_fields );
#endif


    free( EB_WAVE );
    free( J_B0 );
    free( n_e );
    free( data2save );
    return EXIT_SUCCESS;
}//}}}


int make_antenna_profile( gridConfiguration *gridCfg, beamConfiguration *beamCfg, 
                          int exc_signal, 
                          double z2waist,
                          double antField_xy[gridCfg->Nx/2][gridCfg->Ny/2], double antPhaseTerms[gridCfg->Nx/2][gridCfg->Ny/2] ) {
//{{{
// like make_antenna_profile_3 but with previously missing optional for z2waist
// i.e. allowing now for converging beams with their waist not in the antenna plane

    int
        ii,jj;
    double
        antBeam_z_x, antBeam_z_y,
        antBeam_r_x, antBeam_r_y,
        antBeam_wx, antBeam_wy,

        antPhase_x, antPhase_y,
        antPhaseCurve_xR, antPhaseCurve_yR,
        antPhaseCurve_x, antPhaseCurve_y,
        antPhaseGouy_x, antPhaseGouy_y;

    for (ii=0 ; ii<(gridCfg->Nx/2) ; ++ii) {
        // beam coordinate system
        antBeam_r_x  = ((double)ii-(double)beamCfg->ant_x/2.) * cos(beamCfg->antAngle_zx/180.*M_PI);
        antBeam_z_x  = ((double)ii-(double)beamCfg->ant_x/2.) * sin(beamCfg->antAngle_zx/180.*M_PI) * cos(beamCfg->antAngle_zy/180.*M_PI) + z2waist/2;

        // account for tilted Gauss beam
        // w(z)=w0*sqrt(1+(lambda*z/pi*w0^2)^2)
        antBeam_wx  = beamCfg->ant_w0x*(gridCfg->period/2.) * sqrt( 1. + pow( (gridCfg->period/2)*antBeam_z_x/( M_PI*pow(beamCfg->ant_w0x*(gridCfg->period/2), 2) ) , 2)  );

        // phase variation along beam in atenna plane
        antPhase_x  = antBeam_z_x * 2.*M_PI/(gridCfg->period/2.);

        // phase variation due to curvature of phase fronts
        // radius of curvature of phasefronts: R(z)=z+1/z*(pi*w0^2/lambda)^2
        antPhaseCurve_xR    = antBeam_z_x + 1./(antBeam_z_x + 1e-5) 
                                           *pow( M_PI * pow(beamCfg->ant_w0x*gridCfg->period/2., 2) / (gridCfg->period/2) , 2 );
        antPhaseCurve_x     = pow(antBeam_r_x,2) / (2.*antPhaseCurve_xR) * 2.*M_PI/(gridCfg->period/2);

        for (jj=0 ; jj<(gridCfg->Ny/2) ; ++jj) {
            // beam coordinate system
            antBeam_r_y  = ((double)jj-(double)beamCfg->ant_y/2.) * cos(beamCfg->antAngle_zy/180.*M_PI);
            antBeam_z_y  = ((double)jj-(double)beamCfg->ant_y/2.) * sin(beamCfg->antAngle_zy/180.*M_PI) * cos(beamCfg->antAngle_zx/180.*M_PI) + z2waist/2;
        
            // account for tilted Gauss beam
            // w(z)=w0*sqrt(1+(lambda*z/pi*w0^2)^2)
            antBeam_wy  = beamCfg->ant_w0y*(gridCfg->period/2.) * sqrt( 1. + pow( (gridCfg->period/2.)*antBeam_z_y/( M_PI*pow(beamCfg->ant_w0y*(gridCfg->period/2.), 2) ) , 2)  );

            // envelope of antenna field
            antField_xy[ii][jj] = exp( -1.*pow(antBeam_r_x/antBeam_wx, 2) ) 
                                 *exp( -1.*pow(antBeam_r_y/antBeam_wy, 2) );
            // factor: w0/w(z)
            antField_xy[ii][jj] *= beamCfg->ant_w0x*(gridCfg->period/2)/antBeam_wx * beamCfg->ant_w0y*(gridCfg->period/2)/antBeam_wy;

            // phase variation along beam in atenna plane
            antPhase_y          = antBeam_z_y * 2.*M_PI/(gridCfg->period/2.);

            // phase variation due to curvature of phase fronts
            // radius of curvature of phasefronts: R(z)=z+1/z*(pi*w0^2/lambda)^2
            antPhaseCurve_yR    = antBeam_z_y + 1./(antBeam_z_y + 1e-5) 
                                               *pow( M_PI * pow(beamCfg->ant_w0y*gridCfg->period/2., 2) / (gridCfg->period/2.) , 2 );
            antPhaseCurve_y     = pow(antBeam_r_y,2) / (2.*antPhaseCurve_yR) * 2.*M_PI/(gridCfg->period/2.);

            // account for the Gouy-phase
            // phase_Gouy = arctan(z/z_R) 
            // with z_R = pi*w_0^2/lambda the Rayleigh range
            antPhaseGouy_x  = atan( gridCfg->period/2.*antBeam_z_x / (M_PI * pow(beamCfg->ant_w0x*gridCfg->period/2., 2) ) );
            antPhaseGouy_y  = atan( gridCfg->period/2.*antBeam_z_y / (M_PI * pow(beamCfg->ant_w0y*gridCfg->period/2., 2) ) );

                //ant_phase   = .0; <<--- extra phase-term

            antPhaseTerms[ii][jj]   = -antPhase_x 
                                      -antPhase_y
                                      -antPhaseCurve_x
                                      -antPhaseCurve_y
                                      -antPhaseGouy_x
                                      -antPhaseGouy_y;
        }
    }

    return EXIT_SUCCESS;
}//}}}


int make_density_profile( gridConfiguration *gridCfg, 
                          int ne_profile, 
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

    if ( ne_profile == 1 ) {
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
    } else if ( ne_profile == 2 ) {
        // linearly increasing profile with k0Ln as slope
        // n_e(z) = m*z 
        // with m = 2*pi/(k0Ln*lambda) 
        //      z = z-starting_position
        //ne_start_z  = (gridCfg->d_absorb + gridCfg->period)/2;
        ne_start_z  = (gridCfg->d_absorb + 224.155)/2;    // benchmark scenario from STEP project: .15m/l_0*period
        if (ne_start_z%2 != 0)
            ne_start_z  += 1;
        ne_k0Ln     = cntrl_para;
        printf( "make_density_profile: ne_profile = %d, ne_start_z = %ld, k0Ln = %f\n", 
                ne_profile, ne_start_z, ne_k0Ln );
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
    } else if ( ne_profile == 3 ) {
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
    } else if ( ne_profile == 4 ) {
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


int set_densityInAbsorber_v2( double period, int d_absorb, 
                              char absorber[], 
                              size_t N_x, size_t N_y, size_t N_z,
                              double n_e[N_x/2][N_y/2][N_z/2] ) {
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
    ne_dist     = round( period/1 );

    x0          = (double)d_absorb + ne_dist;
    x1          = (double)N_x - (d_absorb + ne_dist);
    y0          = (double)d_absorb + ne_dist;
    y1          = (double)N_y - (d_absorb + ne_dist);
    z0          = (double)d_absorb + ne_dist;
    z1          = (double)N_z - (d_absorb + ne_dist);

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
        for ( x=0; x<(N_x/2) ; ++x ) {
            scale_fact  = +.5*(    tanh(smooth*(x-x0)) + 1);        // x0 boundary
            //printf( "x1: x=%.1f, scale_fact=%f\n", x, scale_fact) ;
            for ( y=0. ; y<(N_y/2) ; ++y )   {
                for (z=0 ; z<(N_z/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }
    // set density in x1 absorber
    if ( strstr(absorber,"x2") ) {
        for ( x=0; x<(N_x/2) ; ++x ) {
            scale_fact  = +.5*(-1.*tanh(smooth*(x-x1)) + 1);       // x1 boundary
            //printf( "x2: x=%.1f, scale_fact=%f\n", x, scale_fact) ;
            for ( y=0. ; y<(N_y/2) ; ++y )   {
                for (z=0 ; z<(N_z/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }

    // set density in y0 absorber
    if ( strstr(absorber,"y1") ) {
        for ( y=0; y<(N_y/2) ; ++y ) {
            scale_fact  = +.5*(    tanh(smooth*(y-y0)) + 1);        // y0 boundary
            //printf( "y1: y=%.1f, scale_fact=%f\n", y, scale_fact) ;
            for ( x=0; x<(N_x/2) ; ++x ) {
                for (z=0 ; z<(N_z/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }
    // set density in y1 absorber
    if ( strstr(absorber,"y2") ) {
        for ( y=0; y<(N_y/2) ; ++y ) {
            scale_fact  = +.5*(-1.*tanh(smooth*(y-y1)) + 1);       // y1 boundary
            //printf( "y2: y=%.1f, scale_fact=%f\n", y, scale_fact) ;
            for ( x=0; x<(N_x/2) ; ++x ) {
                for (z=0 ; z<(N_z/2) ; ++z) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }

    // set density in z0 absorber
    if ( strstr(absorber,"z1") ) {
        for ( z=0 ; z<(N_z/2) ; ++z) {
            scale_fact  = +.5*(    tanh(smooth*(z-z0)) + 1);        // z0 boundary
            //printf( "z1: z=%.1f, scale_fact=%f\n", z, scale_fact) ;
            for ( x=0; x<(N_x/2) ; ++x ) {
                for ( y=0; y<(N_y/2) ; ++y ) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }
    // set density in z1 absorber
    if ( strstr(absorber,"z2") ) {
        for ( z=0 ; z<(N_z/2) ; ++z) {
            scale_fact  = +.5*(-1.*tanh(smooth*(z-z1)) + 1);       // z1 boundary
            //printf( "z2: z=%.1f, scale_fact=%f\n", z, scale_fact) ;
            for ( x=0; x<(N_x/2) ; ++x ) {
                for ( y=0; y<(N_y/2) ; ++y ) {
                    n_e[(int)x][(int)y][(int)z]  *= scale_fact;
                }
            }
        }
    }

    return EXIT_SUCCESS;
} //}}}


int make_B0_profile( gridConfiguration *gridCfg, 
                     int B0_profile, 
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
    if ( B0_profile == 1 ) {
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


int add_source( beamConfiguration *beamCfg, 
                int exc_signal,
                size_t N_x, size_t N_y, size_t N_z,  
                double period, 
                int t_int, double omega_t, 
                double antField_xy[N_x/2][N_y/2], double antPhaseTerms[N_x/2][N_y/2],
                double EB_WAVE[N_x][N_y][N_z] ) {
//{{{

    size_t
        ii, jj;
    double
        t_rise, 
        source;

    if ( exc_signal == 1 ) {
        t_rise  = 1. - exp( -1*pow( ((double)(t_int)/period), 2 )/100. );
#pragma omp parallel for collapse(2) default(shared) private(ii, jj, source)
        for ( ii=2 ; ii<N_x ; ii+=2 ) {
            for ( jj=2 ; jj<N_y ; jj+=2 ) {
                // note: for X-mode injection, switch cos and sin of source_1 and source_2
                //source      = sin(omega_t - aux - curve + GouyPhase_beam + ant_phase/180.*M_PI ) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                // Ex
                EB_WAVE[ii+1][jj  ][beamCfg->ant_z]   += source;
            }
        }
    } else if ( exc_signal == 2) {
        t_rise  = 1. - exp( -1*pow( ((double)(t_int)/period), 2 )/100. );
#pragma omp parallel for collapse(2) default(shared) private(ii, jj, source)
        for ( ii=2 ; ii<N_x ; ii+=2 ) {
            for ( jj=2 ; jj<N_y ; jj+=2 ) {
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                // Bx
                EB_WAVE[ii  ][jj+1][beamCfg->ant_z+1]   += source;
            }
        }
    } else if ( exc_signal == 3) {
        t_rise  = 1. - exp( -1*pow( ((double)(t_int)/period), 2 )/100. );
#pragma omp parallel for collapse(2) default(shared) private(ii, jj, source)
        for ( ii=2 ; ii<N_x ; ii+=2 ) {
            for ( jj=2 ; jj<N_y ; jj+=2 ) {
                // note: for X-mode injection, switch cos and sin of source_1 and source_2
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                // Ex
                EB_WAVE[ii+1][jj  ][beamCfg->ant_z]   += source;
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)] + M_PI/2.) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                // Bx
                EB_WAVE[ii  ][jj+1][beamCfg->ant_z+1] += source*(1.41);
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
    // z2 absorber: z=d_absorb...NZ
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=(gridCfg->Nz-gridCfg->d_absorb) ; kk<gridCfg->Nz-2 ; kk+=2) {      //NZ-d_absorb-2 ???
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
    // x2 absorber: x=d_absorb...NX
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {  
            for (ii=(gridCfg->Nx-gridCfg->d_absorb) ; ii<gridCfg->Nx-2 ; ii+=2) {    //NX-d_absorb-2 ???
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
    // y2 absorber: y=d_absorb...NY
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz-2 ; kk+=2) {
            for (jj=(gridCfg->Ny-gridCfg->d_absorb) ; jj<gridCfg->Ny-2 ; jj+=2) {  //NY-d_absorb-2 ???
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
    // z2 absorber: z=d_absorb...NZ
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
            for (kk=(gridCfg->Nz_ref-gridCfg->d_absorb) ; kk<gridCfg->Nz_ref-2 ; kk+=2) {      //NZ-d_absorb-2 ???
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
    // x2 absorber: x=d_absorb...NX
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (jj=2 ; jj<gridCfg->Ny-2 ; jj+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {  
            for (ii=(gridCfg->Nx-gridCfg->d_absorb) ; ii<gridCfg->Nx-2 ; ii+=2) {    //NX-d_absorb-2 ???
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
    // y2 absorber: y=d_absorb...NY
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
    for (ii=2 ; ii<gridCfg->Nx-2 ; ii+=2) {
        for (kk=2 ; kk<gridCfg->Nz_ref-2 ; kk+=2) {
            for (jj=(gridCfg->Ny-gridCfg->d_absorb) ; jj<gridCfg->Ny-2 ; jj+=2) {  //NY-d_absorb-2 ???
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
    // z2 absorber: z=d_absorb...NZ
    if ( strstr(absorber,"z2") ) {      
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (jj=2 ; jj<N_y-2 ; jj+=2) {
                for (kk=(N_z-d_absorb) ; kk<N_z-2 ; kk+=2) {      //NZ-d_absorb-2 ???
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
    // x2 absorber: x=d_absorb...NX
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (jj=2 ; jj<N_y-2 ; jj+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {  
                for (ii=(N_x-d_absorb) ; ii<N_x-2 ; ii+=2) {    //NX-d_absorb-2 ???
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
    // y2 absorber: y=d_absorb...NY
    if ( strstr(absorber,"y2") ) {
#pragma omp parallel for default(shared) private(ii,jj,kk,damp)
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {
                for (jj=(N_y-d_absorb) ; jj<N_y-2 ; jj+=2) {  //NY-d_absorb-2 ???
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


int abc_Mur_saveOldE_xdir( size_t N_x, size_t N_y, size_t N_z, 
                           double EB_WAVE[N_x][N_y][N_z], 
                           double E_old[8][N_y][N_z] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        jj, kk, 
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj<N_y-2 ; jj+=2) {
        for (kk=2 ; kk<N_z-2 ; kk+=2) {
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
            E_old[4+1][jj  ][kk  ]  = EB_WAVE[N_x-4-offset+1][jj  ][kk  ];
            E_old[6+1][jj  ][kk  ]  = EB_WAVE[N_x-2-offset+1][jj  ][kk  ];
            // Ey: even-odd-even
            E_old[4  ][jj+1][kk  ]  = EB_WAVE[N_x-4-offset  ][jj+1][kk  ];
            E_old[6  ][jj+1][kk  ]  = EB_WAVE[N_x-2-offset  ][jj+1][kk  ];
            // Ez: even-even-odd
            E_old[4  ][jj  ][kk+1]  = EB_WAVE[N_x-4-offset  ][jj  ][kk+1];
            E_old[6  ][jj  ][kk+1]  = EB_WAVE[N_x-2-offset  ][jj  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldE_ydir( size_t N_x, size_t N_y, size_t N_z, 
                           double EB_WAVE[N_x][N_y][N_z], 
                           double E_old[N_y][8][N_z] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, kk,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii<N_x-2 ; ii+=2) {
        for (kk=2 ; kk<N_z-2 ; kk+=2) {
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
            E_old[ii+1][4  ][kk  ]  = EB_WAVE[ii+1][N_y-4-offset  ][kk  ];
            E_old[ii+1][6  ][kk  ]  = EB_WAVE[ii+1][N_y-2-offset  ][kk  ];
            // Ey: even-odd-even
            E_old[ii  ][4+1][kk  ]  = EB_WAVE[ii  ][N_y-4-offset+1][kk  ];
            E_old[ii  ][6+1][kk  ]  = EB_WAVE[ii  ][N_y-2-offset+1][kk  ];
            // Ez: even-even-odd
            E_old[ii  ][4  ][kk+1]  = EB_WAVE[ii  ][N_y-4-offset  ][kk+1];
            E_old[ii  ][6  ][kk+1]  = EB_WAVE[ii  ][N_y-2-offset  ][kk+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abc_Mur_saveOldE_zdir( size_t N_x, size_t N_y, size_t N_z, 
                           double EB_WAVE[N_x][N_y][N_z], 
                           double E_old[N_x][N_y][8] ) {
//{{{

    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj,
        offset;

    offset  = 2;

#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii<N_x-2 ; ii+=2) {
        for (jj=2 ; jj<N_y-2 ; jj+=2) {
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
            E_old[ii+1][jj  ][4  ]  = EB_WAVE[ii+1][jj  ][N_z-4-offset  ];
            E_old[ii+1][jj  ][6  ]  = EB_WAVE[ii+1][jj  ][N_z-2-offset  ];
            // Ey: even-odd-even
            E_old[ii  ][jj+1][4  ]  = EB_WAVE[ii  ][jj+1][N_z-4-offset  ];
            E_old[ii  ][jj+1][6  ]  = EB_WAVE[ii  ][jj+1][N_z-2-offset  ];
            // Ez: even-even-odd
            E_old[ii  ][jj  ][4+1]  = EB_WAVE[ii  ][jj  ][N_z-4-offset+1];
            E_old[ii  ][jj  ][6+1]  = EB_WAVE[ii  ][jj  ][N_z-2-offset+1];
        }
    }
 
    return EXIT_SUCCESS;

}//}}}


int abs_Mur_1st( size_t N_x, size_t N_y, size_t N_z,
                 double dt, double dx, 
                 double EB_WAVE[N_x][N_y][N_z], 
                 double E_old_xdir[8][N_y][N_z], double E_old_ydir[N_x][8][N_z], double E_old_zdir[N_x][N_y][8] ) {
//{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk,
        offset;             // refers to EB_WAVE only

    double
        cnst;

    cnst    = (dt-dx)/(dt+dx);
    offset  = 2;

    // absorber into x-direction
#pragma omp parallel for collapse(2) default(shared) private(jj,kk)
    for (jj=2 ; jj<N_y-2 ; jj+=2) {
        for (kk=2 ; kk<N_z-2 ; kk+=2) {
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
            EB_WAVE[N_x-2-offset+1][jj  ][kk  ]    = E_old_xdir[4+1][jj  ][kk  ]
                + cnst * (    EB_WAVE[N_x-4-offset+1][jj  ][kk  ] 
                          -E_old_xdir[6+1           ][jj  ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[N_x-2-offset  ][jj+1][kk  ]    = E_old_xdir[4  ][jj+1][kk  ]
                + cnst * (    EB_WAVE[N_x-4-offset  ][jj+1][kk  ] 
                          -E_old_xdir[6             ][jj+1][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[N_x-2-offset  ][jj  ][kk+1]    = E_old_xdir[4  ][jj  ][kk+1]
                + cnst * (    EB_WAVE[N_x-4-offset  ][jj  ][kk+1] 
                          -E_old_xdir[6             ][jj  ][kk+1] );
        }
    }

    // absorber into y-direction
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
    for (ii=2 ; ii<N_x-2 ; ii+=2) {
        for (kk=2 ; kk<N_z-2 ; kk+=2) {
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
            EB_WAVE[ii+1][N_y-2-offset  ][kk  ] = E_old_ydir[ii+1][4  ][kk  ]
                + cnst * (    EB_WAVE[ii+1][N_y-4-offset  ][kk  ]
                          -E_old_ydir[ii+1][6      ][kk  ] );
            // Ey: even-odd-even
            EB_WAVE[ii  ][N_y-2-offset+1][kk  ] = E_old_ydir[ii  ][4+1][kk  ]
                + cnst * (    EB_WAVE[ii  ][N_y-4-offset+1][kk  ]
                          -E_old_ydir[ii  ][6+1           ][kk  ] );
            // Ez: even-even-odd
            EB_WAVE[ii  ][N_y-2-offset  ][kk+1] = E_old_ydir[ii  ][4  ][kk+1]
                + cnst * (    EB_WAVE[ii  ][N_y-4-offset  ][kk+1]
                          -E_old_ydir[ii  ][6             ][kk+1] );
        }
    }

    // absorber into z-direction
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
    for (ii=2 ; ii<N_x-2 ; ii+=2) {
        for (jj=2 ; jj<N_y-2 ; jj+=2) {
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
            EB_WAVE[ii+1][jj  ][N_z-2-offset  ]    = E_old_zdir[ii+1][jj  ][4  ]
                + cnst * (    EB_WAVE[ii+1][jj  ][N_z-4-offset  ]
                          -E_old_zdir[ii+1][jj  ][6  ]     );
            // Ey: even-odd-even
            EB_WAVE[ii  ][jj+1][N_z-2-offset  ]    = E_old_zdir[ii  ][jj+1][4  ]
                + cnst * (    EB_WAVE[ii  ][jj+1][N_z-4-offset  ]
                          -E_old_zdir[ii  ][jj+1][6  ]     );
            // Ez: even-even-odd
            EB_WAVE[ii  ][jj  ][N_z-2-offset+1]    = E_old_zdir[ii  ][jj  ][4+1]
                + cnst * (    EB_WAVE[ii  ][jj  ][N_z-4-offset+1]
                          -E_old_zdir[ii  ][jj  ][6+1]     );
        }
    }

    return EXIT_SUCCESS;

} //}}}


int abs_Mur_1st_v2( size_t N_x, size_t N_y, size_t N_z,
                    double dt, double dx, 
                    char absorber[],
                    double EB_WAVE[N_x][N_y][N_z], 
                    double E_old_xdir[8][N_y][N_z], double E_old_ydir[N_x][8][N_z], double E_old_zdir[N_x][N_y][8] ) {
//{{{
    // Ex: odd-even-even
    // Ey: even-odd-even
    // Ez: even-even-odd

    size_t
        ii, jj, kk,
        offset;             // refers to EB_WAVE only

    double
        cnst;

    cnst    = (dt-dx)/(dt+dx);
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
        for (jj=2 ; jj<N_y-2 ; jj+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {
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
        for (jj=2 ; jj<N_y-2 ; jj+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {
                // absorber at x=Nx grid boundary
                // Ex: odd-even-even
                EB_WAVE[N_x-2-offset+1][jj  ][kk  ]    = E_old_xdir[4+1][jj  ][kk  ]
                    + cnst * (    EB_WAVE[N_x-4-offset+1][jj  ][kk  ] 
                              -E_old_xdir[6+1           ][jj  ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[N_x-2-offset  ][jj+1][kk  ]    = E_old_xdir[4  ][jj+1][kk  ]
                    + cnst * (    EB_WAVE[N_x-4-offset  ][jj+1][kk  ] 
                              -E_old_xdir[6             ][jj+1][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[N_x-2-offset  ][jj  ][kk+1]    = E_old_xdir[4  ][jj  ][kk+1]
                    + cnst * (    EB_WAVE[N_x-4-offset  ][jj  ][kk+1] 
                              -E_old_xdir[6             ][jj  ][kk+1] );
            }
        }
    }

    // absorber into y-direction
    if ( strstr(absorber,"y1") ) {
        //printf("abs_Mur_1st_v2: y1\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,kk)
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {
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
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (kk=2 ; kk<N_z-2 ; kk+=2) {
                // absorber at y=Ny grid boundary
                // Ex: odd-even-even
                EB_WAVE[ii+1][N_y-2-offset  ][kk  ] = E_old_ydir[ii+1][4  ][kk  ]
                    + cnst * (    EB_WAVE[ii+1][N_y-4-offset  ][kk  ]
                              -E_old_ydir[ii+1][6      ][kk  ] );
                // Ey: even-odd-even
                EB_WAVE[ii  ][N_y-2-offset+1][kk  ] = E_old_ydir[ii  ][4+1][kk  ]
                    + cnst * (    EB_WAVE[ii  ][N_y-4-offset+1][kk  ]
                              -E_old_ydir[ii  ][6+1           ][kk  ] );
                // Ez: even-even-odd
                EB_WAVE[ii  ][N_y-2-offset  ][kk+1] = E_old_ydir[ii  ][4  ][kk+1]
                    + cnst * (    EB_WAVE[ii  ][N_y-4-offset  ][kk+1]
                              -E_old_ydir[ii  ][6             ][kk+1] );
            }
        }
    }

    // absorber into z-direction
    if ( strstr(absorber,"z1") ) {
        //printf("abs_Mur_1st_v2: z1\n");
#pragma omp parallel for collapse(2) default(shared) private(ii,jj)
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (jj=2 ; jj<N_y-2 ; jj+=2) {
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
        for (ii=2 ; ii<N_x-2 ; ii+=2) {
            for (jj=2 ; jj<N_y-2 ; jj+=2) {
                // absorber at z=Nz grid boundary
                // Ex: odd-even-even
                EB_WAVE[ii+1][jj  ][N_z-2-offset  ]    = E_old_zdir[ii+1][jj  ][4  ]
                    + cnst * (    EB_WAVE[ii+1][jj  ][N_z-4-offset  ]
                              -E_old_zdir[ii+1][jj  ][6  ]     );
                // Ey: even-odd-even
                EB_WAVE[ii  ][jj+1][N_z-2-offset  ]    = E_old_zdir[ii  ][jj+1][4  ]
                    + cnst * (    EB_WAVE[ii  ][jj+1][N_z-4-offset  ]
                              -E_old_zdir[ii  ][jj+1][6  ]     );
                // Ez: even-even-odd
                EB_WAVE[ii  ][jj  ][N_z-2-offset+1]    = E_old_zdir[ii  ][jj  ][4+1]
                    + cnst * (    EB_WAVE[ii  ][jj  ][N_z-4-offset+1]
                              -E_old_zdir[ii  ][jj  ][6+1]     );
            }
        }
    }

    return EXIT_SUCCESS;

} //}}}


int advance_J( size_t dim1, size_t dim2, size_t dim3, 
               double EB_WAVE[dim1][dim2][dim3], double J_B0[dim1][dim2][dim3],
               double n_e[dim1/2][dim2/2][dim3/2], 
               double dt ) {
//{{{
    size_t
        ii, jj, kk;

    // Jx: odd-even-even
    // Jy: even-odd-even
    // Jz: even-even-odd
    // B0x: even-odd-odd
    // B0y: odd-even-odd
    // B0z: odd-odd-even
//#pragma omp parallel for collapse(2) default(shared) private(k,j, Jx_tmp1,Jy_tmp1,Jz_tmp1, Jx_tmp2,Jy_tmp2,Jz_tmp2, Jx_tmp3,Jy_tmp3,Jz_tmp3 ) // collapse ???
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii<dim1-2 ; ii+=2) {
        for (jj=2 ; jj<dim2-2 ; jj+=2) {
            for (kk=2 ; kk<dim3-2 ; kk+=2) {
                // Jx: odd-even-even
                J_B0[ii+1][jj  ][kk  ]    += + dt*(
                        pow(2*M_PI,2) * n_e[(ii/2)][(jj/2)][(kk/2)] * EB_WAVE[ii+1][jj  ][kk  ]
                        - 2*M_PI * ( J_B0[ii  ][jj+1][kk  ]*J_B0[ii+1][jj+1][kk  ]        // +Jy*B0z
                                    -J_B0[ii  ][jj  ][kk+1]*J_B0[ii+1][jj  ][kk+1]        // -Jz*B0y
                              )
                        );
                // Jy: even-odd-even
                J_B0[ii  ][jj+1][kk  ]    += + dt*(
                        pow(2*M_PI,2) * n_e[(ii/2)][(jj/2)][(kk/2)] * EB_WAVE[ii  ][jj+1][kk  ]
                        -2*M_PI * (-J_B0[ii+1][jj  ][kk  ]*J_B0[ii+1][jj+1][kk  ]         // -Jx*B0z
                                   +J_B0[ii  ][jj  ][kk+1]*J_B0[ii  ][jj+1][kk+1]         // +Jz*B0x
                              )
                        );
                // Jz: even-even-odd
                J_B0[ii  ][jj  ][kk+1]    += + dt*(
                        pow(2*M_PI,2) * n_e[(ii/2)][(jj/2)][(kk/2)] * EB_WAVE[ii  ][jj  ][kk+1]
                        -2*M_PI * ( J_B0[ii+1][jj  ][kk  ]*J_B0[ii+1][jj  ][kk+1]        // +Jx*B0y
                                   -J_B0[ii  ][jj+1][kk  ]*J_B0[ii  ][jj+1][kk+1]        // -Jy*B0x
                              )
                        );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_B( size_t dim1, size_t dim2, size_t dim3, 
               double EB_WAVE[dim1][dim2][dim3], 
               double dx, double dt ) {
//{{{
    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk) 
    for (ii=2 ; ii<dim1-2 ; ii+=2) {
        for (jj=2 ; jj<dim2-2 ; jj+=2) {
            for (kk=2 ; kk<dim3-2 ; kk+=2) {
                // -dBx/dt = dEz/dy - dEy/dz
                EB_WAVE[ii  ][jj+1][kk+1]   += -1.*dt/dx*(
                        +EB_WAVE[ii  ][jj+2][kk+1] - EB_WAVE[ii  ][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+2] + EB_WAVE[ii  ][jj+1][kk  ]
                        );
                // -dBy/dt = dEx/dz - dEz/dx
                EB_WAVE[ii+1][jj  ][kk+1] += -1.*dt/dx*(
                        +EB_WAVE[ii+1][jj  ][kk+2] - EB_WAVE[ii+1][jj  ][kk  ]
                        -EB_WAVE[ii+2][jj  ][kk+1] + EB_WAVE[ii  ][jj  ][kk+1]
                        );
                // -dBz/dt = dEy/dx - dEx/dy
                EB_WAVE[ii+1][jj+1][kk  ] += -1.*dt/dx*(
                        +EB_WAVE[ii+2][jj+1][kk  ] - EB_WAVE[ii  ][jj+1][kk  ]
                        -EB_WAVE[ii+1][jj+2][kk  ] + EB_WAVE[ii+1][jj  ][kk  ]
                        );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_E( size_t dim1, size_t dim2, size_t dim3, 
               double EB_WAVE[dim1][dim2][dim3], double J_B0[dim1][dim2][dim3],
               double dx, double dt ) {
//{{{
    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii<dim1-2 ; ii+=2) {
        for (jj=2 ; jj<dim2-2 ; jj+=2) {
            for (kk=2 ; kk<dim3-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                EB_WAVE[ii+1][jj  ][kk  ] += dt/dx*(
                        +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                        -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                        ) - dt*J_B0[ii+1][jj  ][kk  ];
                // dEy/dt = (dBx/dz - dBz/dx)
                EB_WAVE[ii  ][jj+1][kk  ] += dt/dx*(
                        +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                        -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                        ) - dt*J_B0[ii  ][jj+1][kk  ];
                // dEz/dt = (dBy/dx - dBx/dy)
                EB_WAVE[ii  ][jj  ][kk+1] += dt/dx*(
                        +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                        ) - dt*J_B0[ii  ][jj  ][kk+1];
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int advance_E_vacuum( size_t dim1, size_t dim2, size_t dim3, 
                      double EB_WAVE[dim1][dim2][dim3],
                      double dx, double dt ) {
//{{{
    size_t
        ii, jj, kk;

#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=2 ; ii<dim1-2 ; ii+=2) {
        for (jj=2 ; jj<dim2-2 ; jj+=2) {
            for (kk=2 ; kk<dim3-2 ; kk+=2) {
                // dEx/dt = (dBz/dy - dBy/dz)
                EB_WAVE[ii+1][jj  ][kk  ] += dt/dx*(
                        +EB_WAVE[ii+1][jj+1][kk  ] - EB_WAVE[ii+1][jj-1][kk  ]
                        -EB_WAVE[ii+1][jj  ][kk+1] + EB_WAVE[ii+1][jj  ][kk-1]
                        );
                // dEy/dt = (dBx/dz - dBz/dx)
                EB_WAVE[ii  ][jj+1][kk  ] += dt/dx*(
                        +EB_WAVE[ii  ][jj+1][kk+1] - EB_WAVE[ii  ][jj+1][kk-1]
                        -EB_WAVE[ii+1][jj+1][kk  ] + EB_WAVE[ii-1][jj+1][kk  ]
                        );
                // dEz/dt = (dBy/dx - dBx/dy)
                EB_WAVE[ii  ][jj  ][kk+1] += dt/dx*(
                        +EB_WAVE[ii+1][jj  ][kk+1] - EB_WAVE[ii-1][jj  ][kk+1]
                        -EB_WAVE[ii  ][jj+1][kk+1] + EB_WAVE[ii  ][jj-1][kk+1]
                        );
            }
        }
    }
    return EXIT_SUCCESS;
}//}}}


int calc_poynt_1( size_t N_x, size_t N_y, size_t N_z, 
                  int pwr_dect, char absorber[], double poynt[3],  
                  double EB_WAVE[N_x][N_y][N_z] ) {
//{{{

    size_t
        ii, jj, kk;
    double
        poynt_x, poynt_y, poynt_z;

    poynt_x = .0;
    poynt_y = .0;
    poynt_z = .0;

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
    
    if ( strcmp(absorber,"z1") == 0 ) {
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                poynt_x += ( EB_WAVE[ii  ][jj+1][pwr_dect  ]
                            *EB_WAVE[ii+1][jj+1][pwr_dect  ]
                            -EB_WAVE[ii  ][jj  ][pwr_dect+1]
                            *EB_WAVE[ii+1][jj  ][pwr_dect+1] );
                poynt_y += ( EB_WAVE[ii  ][jj  ][pwr_dect+1]
                            *EB_WAVE[ii  ][jj+1][pwr_dect+1]
                            -EB_WAVE[ii+1][jj  ][pwr_dect  ]
                            *EB_WAVE[ii+1][jj+1][pwr_dect  ] );
                poynt_z += ( EB_WAVE[ii+1][jj  ][pwr_dect  ]
                            *EB_WAVE[ii+1][jj  ][pwr_dect+1]
                            -EB_WAVE[ii  ][jj+1][pwr_dect  ]
                            *EB_WAVE[ii  ][jj+1][pwr_dect+1] );
            }
        }
    } else if ( strcmp(absorber,"z2") == 0 ) {
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z2-plane
                poynt_x += ( EB_WAVE[ii  ][jj+1][N_z-pwr_dect  ]
                            *EB_WAVE[ii+1][jj+1][N_z-pwr_dect  ]
                            -EB_WAVE[ii  ][jj  ][N_z-pwr_dect-1]
                            *EB_WAVE[ii+1][jj  ][N_z-pwr_dect-1] );
                poynt_y += ( EB_WAVE[ii  ][jj  ][N_z-pwr_dect-1]
                            *EB_WAVE[ii  ][jj+1][N_z-pwr_dect-1]
                            -EB_WAVE[ii+1][jj  ][N_z-pwr_dect  ]
                            *EB_WAVE[ii+1][jj+1][N_z-pwr_dect  ] );
                poynt_z += ( EB_WAVE[ii+1][jj  ][N_z-pwr_dect  ]
                            *EB_WAVE[ii+1][jj  ][N_z-pwr_dect-1]
                            -EB_WAVE[ii  ][jj+1][N_z-pwr_dect  ]
                            *EB_WAVE[ii  ][jj+1][N_z-pwr_dect-1] );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x1-plane
                poynt_x += ( EB_WAVE[pwr_dect  ][jj+1][kk  ]
                            *EB_WAVE[pwr_dect+1][jj+1][kk  ]
                            -EB_WAVE[pwr_dect  ][jj  ][kk+1]
                            *EB_WAVE[pwr_dect+1][jj  ][kk+1] );
                poynt_y += ( EB_WAVE[pwr_dect  ][jj  ][kk+1]
                            *EB_WAVE[pwr_dect  ][jj+1][kk+1]
                            -EB_WAVE[pwr_dect+1][jj  ][kk  ]
                            *EB_WAVE[pwr_dect+1][jj+1][kk  ] );
                poynt_z += ( EB_WAVE[pwr_dect+1][jj  ][kk  ]
                            *EB_WAVE[pwr_dect+1][jj  ][kk+1]
                            -EB_WAVE[pwr_dect  ][jj+1][kk  ]
                            *EB_WAVE[pwr_dect  ][jj+1][kk+1] );
            }
        }
    } else if ( strcmp(absorber,"x2") == 0 ) {
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x2-plane
                poynt_x += ( EB_WAVE[N_x-pwr_dect  ][jj+1][kk  ]
                            *EB_WAVE[N_x-pwr_dect-1][jj+1][kk  ]
                            -EB_WAVE[N_x-pwr_dect  ][jj  ][kk+1]
                            *EB_WAVE[N_x-pwr_dect-1][jj  ][kk+1] );
                poynt_y += ( EB_WAVE[N_x-pwr_dect  ][jj  ][kk+1]
                            *EB_WAVE[N_x-pwr_dect  ][jj+1][kk+1]
                            -EB_WAVE[N_x-pwr_dect-1][jj  ][kk  ]
                            *EB_WAVE[N_x-pwr_dect-1][jj+1][kk  ] );
                poynt_z += ( EB_WAVE[N_x-pwr_dect-1][jj  ][kk  ]
                            *EB_WAVE[N_x-pwr_dect-1][jj  ][kk+1]
                            -EB_WAVE[N_x-pwr_dect  ][jj+1][kk  ]
                            *EB_WAVE[N_x-pwr_dect  ][jj+1][kk+1] );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y1-plane
                poynt_x += ( EB_WAVE[ii  ][pwr_dect+1][kk  ]
                            *EB_WAVE[ii+1][pwr_dect+1][kk  ]
                            -EB_WAVE[ii  ][pwr_dect  ][kk+1]
                            *EB_WAVE[ii+1][pwr_dect  ][kk+1] );
                poynt_y += ( EB_WAVE[ii  ][pwr_dect  ][kk+1]
                            *EB_WAVE[ii  ][pwr_dect+1][kk+1]
                            -EB_WAVE[ii+1][pwr_dect  ][kk  ]
                            *EB_WAVE[ii+1][pwr_dect+1][kk  ] );
                poynt_z += ( EB_WAVE[ii+1][pwr_dect  ][kk  ]
                            *EB_WAVE[ii+1][pwr_dect  ][kk+1]
                            -EB_WAVE[ii  ][pwr_dect+1][kk  ]
                            *EB_WAVE[ii  ][pwr_dect+1][kk+1] );
            }
        }
    } else if ( strcmp(absorber,"y2") == 0 ) {
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y2-plane
                poynt_x += ( EB_WAVE[ii  ][N_y-pwr_dect-1][kk  ]
                            *EB_WAVE[ii+1][N_y-pwr_dect-1][kk  ]
                            -EB_WAVE[ii  ][N_y-pwr_dect  ][kk+1]
                            *EB_WAVE[ii+1][N_y-pwr_dect  ][kk+1] );
                poynt_y += ( EB_WAVE[ii  ][N_y-pwr_dect  ][kk+1]
                            *EB_WAVE[ii  ][N_y-pwr_dect-1][kk+1]
                            -EB_WAVE[ii+1][N_y-pwr_dect  ][kk  ]
                            *EB_WAVE[ii+1][N_y-pwr_dect-1][kk  ] );
                poynt_z += ( EB_WAVE[ii+1][N_y-pwr_dect  ][kk  ]
                            *EB_WAVE[ii+1][N_y-pwr_dect  ][kk+1]
                            -EB_WAVE[ii  ][N_y-pwr_dect-1][kk  ]
                            *EB_WAVE[ii  ][N_y-pwr_dect-1][kk+1] );
            }
        }
    }

    poynt[0]    = poynt_x;
    poynt[1]    = poynt_y;
    poynt[2]    = poynt_z;

    return EXIT_SUCCESS;
}//}}}


double calc_poynt_2( size_t N_x, size_t N_y, size_t N_z, 
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z] ) {
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
    
    if ( strcmp(absorber,"z1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE[ii+1][jj  ][pwr_dect  ]
                          *EB_WAVE[ii+1][jj  ][pwr_dect+1]
                          -EB_WAVE[ii  ][jj+1][pwr_dect  ]
                          *EB_WAVE[ii  ][jj+1][pwr_dect+1] );
            }
        }
    } else if ( strcmp(absorber,"z2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z2-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE[ii+1][jj  ][N_z-pwr_dect  ]
                          *EB_WAVE[ii+1][jj  ][N_z-pwr_dect-1]
                          -EB_WAVE[ii  ][jj+1][N_z-pwr_dect  ]
                          *EB_WAVE[ii  ][jj+1][N_z-pwr_dect-1] );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
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
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x2-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[N_x-pwr_dect  ][jj+1][kk  ]
                          *EB_WAVE[N_x-pwr_dect-1][jj+1][kk  ]
                          -EB_WAVE[N_x-pwr_dect  ][jj  ][kk+1]
                          *EB_WAVE[N_x-pwr_dect-1][jj  ][kk+1] );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
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
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y2-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][N_y-pwr_dect  ][kk+1]
                          *EB_WAVE[ii  ][N_y-pwr_dect-1][kk+1]
                          -EB_WAVE[ii+1][N_y-pwr_dect  ][kk  ]
                          *EB_WAVE[ii+1][N_y-pwr_dect-1][kk  ] );
            }
        }
    }
    
    return poynt;
}//}}}


double calc_poynt_3( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] ) {
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
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
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
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
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
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z2-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE[ii+1][jj  ][N_z-pwr_dect  ]
                          *EB_WAVE[ii+1][jj  ][N_z-pwr_dect-1]
                          -EB_WAVE[ii  ][jj+1][N_z-pwr_dect  ]
                          *EB_WAVE[ii  ][jj+1][N_z-pwr_dect-1] );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
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
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x2-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[N_x-pwr_dect  ][jj+1][kk  ]
                          *EB_WAVE[N_x-pwr_dect-1][jj+1][kk  ]
                          -EB_WAVE[N_x-pwr_dect  ][jj  ][kk+1]
                          *EB_WAVE[N_x-pwr_dect-1][jj  ][kk+1] );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
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
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y2-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][N_y-pwr_dect  ][kk+1]
                          *EB_WAVE[ii  ][N_y-pwr_dect-1][kk+1]
                          -EB_WAVE[ii+1][N_y-pwr_dect  ][kk  ]
                          *EB_WAVE[ii+1][N_y-pwr_dect-1][kk  ] );
            }
        }
    }
    
    return poynt;
} //}}}


double calc_poynt_4( size_t N_x, size_t N_y, size_t N_z, size_t N_z_ref,
                     int pwr_dect, char absorber[],
                     double EB_WAVE[N_x][N_y][N_z], double EB_WAVE_ref[N_x][N_y][N_z_ref] ) {
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
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
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
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z1-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( ( EB_WAVE[ii+1][jj  ][pwr_dect  ] - EB_WAVE_ref[ii+1][jj  ][pwr_dect  ] )
                          *( EB_WAVE[ii+1][jj  ][pwr_dect+1] - EB_WAVE_ref[ii+1][jj  ][pwr_dect+1] )
                          -( EB_WAVE[ii  ][jj+1][pwr_dect  ] - EB_WAVE_ref[ii  ][jj+1][pwr_dect  ] ) 
                          *( EB_WAVE[ii  ][jj+1][pwr_dect+1] - EB_WAVE_ref[ii  ][jj+1][pwr_dect+1] ) );
                //poynt += ( ( -EB_WAVE[ii+1][jj  ][pwr_dect  ] + EB_WAVE_ref[ii+1][jj  ][pwr_dect  ] )
                //          *( -EB_WAVE[ii+1][jj  ][pwr_dect+1] + EB_WAVE_ref[ii+1][jj  ][pwr_dect+1] )
                //          -( -EB_WAVE[ii  ][jj+1][pwr_dect  ] + EB_WAVE_ref[ii  ][jj+1][pwr_dect  ] ) 
                //          *( -EB_WAVE[ii  ][jj+1][pwr_dect+1] + EB_WAVE_ref[ii  ][jj+1][pwr_dect+1] ) );
            }
        }
    } else if ( strcmp(absorber,"z2") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,jj) reduction(+:poynt)
        for (ii=pwr_dect ; ii<(N_x-pwr_dect-2) ; ii+=2) {
            for (jj=pwr_dect ; jj<(N_y-pwr_dect-2) ; jj+=2) {
                // z2-plane
                // Pz = Ex*Hy - Ey*Hx
                poynt += ( EB_WAVE[ii+1][jj  ][N_z-pwr_dect  ]
                          *EB_WAVE[ii+1][jj  ][N_z-pwr_dect-1]
                          -EB_WAVE[ii  ][jj+1][N_z-pwr_dect  ]
                          *EB_WAVE[ii  ][jj+1][N_z-pwr_dect-1] );
                }
            }
    } else if ( strcmp(absorber,"x1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(jj,kk) reduction(+:poynt)
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
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
        for (jj=pwr_dect ; jj<=(N_y-pwr_dect-2) ; jj+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // x2-plane
                // Px = Ey*Hz - Ez*Hy
                poynt += ( EB_WAVE[N_x-pwr_dect  ][jj+1][kk  ]
                          *EB_WAVE[N_x-pwr_dect-1][jj+1][kk  ]
                          -EB_WAVE[N_x-pwr_dect  ][jj  ][kk+1]
                          *EB_WAVE[N_x-pwr_dect-1][jj  ][kk+1] );
                }
            }
    } else if ( strcmp(absorber,"y1") == 0 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii,kk) reduction(+:poynt)
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
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
        for (ii=pwr_dect ; ii<=(N_x-pwr_dect-2) ; ii+=2) {
            for (kk=pwr_dect ; kk<=(N_z-pwr_dect-2) ; kk+=2) {
                // y2-plane
                // Py = Ez*Hx - Ex*Hz
                poynt += ( EB_WAVE[ii  ][N_y-pwr_dect  ][kk+1]
                          *EB_WAVE[ii  ][N_y-pwr_dect-1][kk+1]
                          -EB_WAVE[ii+1][N_y-pwr_dect  ][kk  ]
                          *EB_WAVE[ii+1][N_y-pwr_dect-1][kk  ] );
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
int detAnt1D_storeValues( size_t N_x, size_t N_y, size_t N_z,
                          size_t detAnt_ypos, size_t detAnt_zpos,
                          int tt, double period,  
                          double EB_WAVE[N_x][N_y][N_z], double detAnt_fields[N_x/2][5] ) { 
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
    for ( ii=2 ; ii <= N_x-2 ; ii+=2 ) {
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
        detAnt_fields[ii/2][4]  += ( foo * sqrt(1./( (double)(tt)/(double)(period) + 1e-6 )) );

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
    }

    // apply shuffle and gzip filters, if available
    if (filter_avail) {
        // create dataset creation property list
        dcpl = H5Pcreate( H5P_DATASET_CREATE );
        // add shuffle filter and gzip compression filter
        // note that the order of filter is significant: first shuffle!
        // order of filters applied correspond to order in which they are invoked when writin gdata
        status = H5Pset_shuffle( dcpl );
        status = H5Pset_deflate( dcpl, 9 );
        status = H5Pset_chunk(dcpl, 1, chunk );
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

    status       = H5Dclose(dataset_id_i);
    status       = H5Sclose(dataspace_id_i);


    // store position
    dims[0] = 1;
    data2save[0] = (double)(detAnt_ypos);
    data2save[1] = (double)(detAnt_zpos);
    dataspace_id = H5Screate_simple( 1, dims, NULL); 
    // detAnt_ypos
    printf( "start to create dataset 'detAnt_ypos'\n" );
    dataset_id   = H5Dcreate( group_id__detAnt, "detAnt_ypos", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
    status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data2save[0]);
    status       = H5Dclose(dataset_id);
//    status       = H5Sclose(dataspace_id);
    // detAnt_zpos
//    dataspace_id = H5Screate_simple( 1, dims, NULL); 
    printf( "start to create dataset 'detAnt_zpos'\n" );
    dataset_id   = H5Dcreate( group_id__detAnt, "detAnt_zpos", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
    status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data2save[1]);
    status       = H5Dclose(dataset_id);
    status       = H5Sclose(dataspace_id);

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
    status       = H5Dclose(dataset_id);
//    status       = H5Sclose(dataspace_id);
    // store sum_ExEx

    // store sum_EyEy
    set2zero_1D( N_x/2, data2save );
    for ( ii=2 ; ii<=N_x-2 ; ii+=2 )
        data2save[ii/2] = detAnt_fields[ii/2][1];
//    dataspace_id = H5Screate_simple( 1, dims, NULL);   
    printf( "start to create dataset 'sum_EyEy'\n" );
    if (filter_avail)
        dataset_id   = H5Dcreate( group_id__detAnt, "sum_EyEy", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, dcpl, H5P_DEFAULT);  
    else
        dataset_id   = H5Dcreate( group_id__detAnt, "sum_EyEy", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
    status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2save);
    status       = H5Dclose(dataset_id);
//    status       = H5Sclose(dataspace_id);

    // store sum_EzEz
    set2zero_1D( N_x/2, data2save );
    for ( ii=2 ; ii<=N_x-2 ; ii+=2 )
        data2save[ii/2] = detAnt_fields[ii/2][2];
//    dataspace_id = H5Screate_simple( 1, dims, NULL);   
    printf( "start to create dataset 'sum_EzEz'\n" );
    if (filter_avail)
        dataset_id   = H5Dcreate( group_id__detAnt, "sum_EzEz", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, dcpl, H5P_DEFAULT);  
    else
        dataset_id   = H5Dcreate( group_id__detAnt, "sum_EzEz", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
    status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2save);
    status       = H5Dclose(dataset_id);
//    status       = H5Sclose(dataspace_id);
    
    // store sum_EE
    set2zero_1D( N_x/2, data2save );
    for ( ii=2 ; ii<=N_x-2 ; ii+=2 )
        data2save[ii/2] = detAnt_fields[ii/2][3];
//    dataspace_id = H5Screate_simple( 1, dims, NULL);   
    printf( "start to create dataset 'sum_EE'\n" );
    if (filter_avail)
        dataset_id   = H5Dcreate( group_id__detAnt, "sum_EE", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, dcpl, H5P_DEFAULT);  
    else
        dataset_id   = H5Dcreate( group_id__detAnt, "sum_EE", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
    status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2save);
    status       = H5Dclose(dataset_id);
//    status       = H5Sclose(dataspace_id);

    // store rmsE
    set2zero_1D( N_x/2, data2save );
    for ( ii=2 ; ii<=N_x-2 ; ii+=2 )
        data2save[ii/2] = detAnt_fields[ii/2][4];
//    dataspace_id = H5Screate_simple( 1, dims, NULL);   
    printf( "start to create dataset 'rms_E'\n" );
    if (filter_avail)
        dataset_id   = H5Dcreate( group_id__detAnt, "rms_E", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, dcpl, H5P_DEFAULT);  
    else
        dataset_id   = H5Dcreate( group_id__detAnt, "rms_E", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
    status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2save);
    status       = H5Dclose(dataset_id);

    status       = H5Sclose(dataspace_id);
    
    // terminate access and free ressources/identifiers
    if (filter_avail)
        status = H5Pclose( dcpl );
    status = H5Gclose( group_id__detAnt );
    // file 
    status = H5Fclose(file_id);

    return EXIT_SUCCESS;
}//#}}}
#endif


#ifdef HDF5
int writeMyHDF_v4( int dim0, int dim1, int dim2, char filename[], char dataset[], double array_3D[dim0][dim1][dim2] ) {
    //#{{{

    // hdf related variables
    hid_t       file_id, dataset_id, dataspace_id;      // object identifiers
    hsize_t     dims[3];                                // size used for dimensions
    herr_t      status;                                 // function return value

    // hdf5 related variables used for applying shuffle and compression filter
    hid_t       dcpl;
    hsize_t     chunk[3];
    unsigned int    filter_info;
    int         filter_avail;

    // required for check if hdf5-file already exists
    struct      stat st;

    // assume as default setting, that filters are available
    filter_avail = 1;

    // check if specified hdf5 file already exists
    // if not, create new one; if yes, open and add dataset to it
    if ( stat( filename, &st )==0 ) {
        // open file for read + write access
        file_id = H5Fopen( filename,            // filename
                           H5F_ACC_RDWR,        // allow read & write access (_RDONLY for read only)
                           H5P_DEFAULT);        // file access property list (default one)

        // hdf5 version 1.8.0 introduced H5Lexists to check if link (to group or dataset) exists in hdf5-file
#if H5_VERS_MAJOR>=1 && H5_VERS_MINOR>=8
        if ( H5_VERS_MINOR >= 10 ) {
            printf( "WARNING: hdf5 version 1.10 (or larger is used)\n" );
            printf( "         behavior of H5Lexists was slightly changed in this version\n" );
            printf( "         for details, see https://support.hdfgroup.org/HDF5/doc/RM/RM_H5L.html#Link-Exists\n" );
        }
        if ( H5Lexists( file_id,                // file or group identifier
                        dataset,                // name of link (to group or dataset) to check
                        H5P_DEFAULT )           // link access property list identifiert
                > 0 ) {                         // NOTE: for version 1.8.10, this might be slightly different
            printf( "ERROR: dataset named '%s' already exists in file '%s'\n", dataset, filename );
            printf( "       dataset will NOT be saved (no overwrite by default)\n" );
            status = H5Fclose(file_id);
            return EXIT_FAILURE;
        }
#endif
    } else {
        // create a new file using default properties.
        file_id = H5Fcreate( filename,          // filename
                             H5F_ACC_TRUNC,     // how file should be created (removes existing file)
                             H5P_DEFAULT,       // file creating property list
                             H5P_DEFAULT);      // file access property list
    }

    // create simple data space for the dataset
    // (simple = regular N-dimensional array, i.e. data on regular rectangular grid)
    // (complex = e.g.: irregular grids resulting from dynamic refinement of mesh)
    dims[0] = dim0;
    dims[1] = dim1;
    dims[2] = dim2;
    dataspace_id = H5Screate_simple( 3,     // number of dimensions of dataspace 
                                     dims,  // size of array in each dimension
                                     NULL); // allow stretching of data space (NULL=no)

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
    }

    // apply shuffle and gzip filters, if available
    if (filter_avail) {
        // set chunk size to be same as dimension (might not be optimal, but seems same as h5repack)
        chunk[0] = dim0;
        chunk[1] = dim1;
        chunk[2] = dim2;
        // create dataset creation property list
        dcpl = H5Pcreate( H5P_DATASET_CREATE );
        // add shuffle filter and gzip compression filter
        // note that the order of filter is significant: first shuffle!
        // order of filters applied correspond to order in which they are invoked when writin gdata
        status = H5Pset_shuffle( dcpl );
        status = H5Pset_deflate( dcpl, 9 );
        status = H5Pset_chunk(dcpl,             // dataset creation property list identifier
                              3,                // number of dimensions of each chunk
                              chunk );          // array defining size, in dataset elements, of each chunk
        // create the dataset
        dataset_id = H5Dcreate( file_id,        // file identifier (or group identifier)
                                dataset,        // name of dataset (relative to group specified, if speficied)
                                H5T_NATIVE_DOUBLE, // datatype to use when creating dataset
                                dataspace_id,   // dataspace identifier
                                H5P_DEFAULT,    // link creation property list (was dataset creating property list <=v1.6)
                                dcpl,           // dataset creation property list (added in HDF5v1.8)
                                H5P_DEFAULT);   // dataset access property list (added in HDF5v1.8)
    } else {
        // create the dataset
        dataset_id = H5Dcreate( file_id,        // file identifier (or group identifier)
                                dataset,        // name of dataset (relative to group specified, if speficied)
                                H5T_NATIVE_DOUBLE, // datatype to use when creating dataset
                                dataspace_id,   // dataspace identifier
                                H5P_DEFAULT,    // link creation property list (was dataset creating property list <=v1.6)
                                H5P_DEFAULT,    // dataset creation property list (added in HDF5v1.8)
                                H5P_DEFAULT);   // dataset access property list (added in HDF5v1.8)
    }

    // write the dataset
    status = H5Dwrite( dataset_id,          // dataset identifier
                       H5T_NATIVE_DOUBLE,    // informs hdf about format of data in memory of computer
                       H5S_ALL,             // identifier of memory dataspace
                       H5S_ALL,             // file space identifier
                       H5P_DEFAULT,         // data transfer property list
//                       array_2D[0]);        // pointer to data array
                       array_3D);        // pointer to data array

    // terminate access and free ressources/identifiers
    // dataset creation property list
    if (filter_avail)
        status = H5Pclose(dcpl);
    // dataset
    status = H5Dclose(dataset_id);
    // data space
    status = H5Sclose(dataspace_id);
    // file 
    status = H5Fclose(file_id);
    
    return EXIT_SUCCESS;
}//#}}}
#endif


#ifdef HDF5
int writeConfig2HDF( char filename[], int N_x, int N_y, int N_z, int period, int d_absorb ) {
    //#{{{

    long        data2write_long[1];
    double      data2write_double[1];

    // hdf related variables
    hid_t       file_id, dataset_id, dataspace_id;      // object identifiers
    hsize_t     dims[1];                                // size used for dimensions
    herr_t      status;                                 // function return value

    // note that shuffle and compression filter is not applied here, as only single values are saved

    // required for check if hdf5-file already exists
    struct      stat st;

    // check if specified hdf5 file already exists
    // if not, create new one; if yes, open and add dataset to it
    if ( stat( filename, &st )==0 ) {
        // open file for read + write access
        file_id = H5Fopen( filename,            // filename
                           H5F_ACC_RDWR,        // allow read & write access (_RDONLY for read only)
                           H5P_DEFAULT);        // file access property list (default one)

        // hdf5 version 1.8.0 introduced H5Lexists to check if link (to group or dataset) exists in hdf5-file
#if H5_VERS_MAJOR>=1 && H5_VERS_MINOR>=8
        if ( H5_VERS_MINOR >= 10 ) {
            printf( "WARNING: hdf5 version 1.10 (or larger is used)\n" );
            printf( "         behavior of H5Lexists was slightly changed in this version\n" );
            printf( "         for details, see https://support.hdfgroup.org/HDF5/doc/RM/RM_H5L.html#Link-Exists\n" );
        }
        if ( H5Lexists( file_id,                // file or group identifier
                        "/config",              // name of link (to group or dataset) to check
                        H5P_DEFAULT )           // link access property list identifiert
                > 0 ) {                         // NOTE: for version 1.8.10, this might be slightly different
            printf( "ERROR: dataset named '/config' already exists in file '%s'\n", filename );
            printf( "       dataset will NOT be saved (no overwrite by default)\n" );
            status = H5Fclose(file_id);
            return EXIT_FAILURE;
        }
#endif
    } else {
        // create a new file using default properties.
        file_id = H5Fcreate( filename,          // filename
                             H5F_ACC_TRUNC,     // how file should be created (removes existing file)
                             H5P_DEFAULT,       // file creating property list
                             H5P_DEFAULT);      // file access property list
    }

    // create new group for config parameters
    H5Gcreate( file_id, 
               "/config",
               H5P_DEFAULT,
               H5P_DEFAULT,
               H5P_DEFAULT);

    // create simple data space for the dataset
    // (simple = regular N-dimensional array, i.e. data on regular rectangular grid)
    // (complex = e.g.: irregular grids resulting from dynamic refinement of mesh)
    dims[0] = 1;
    dataspace_id = H5Screate_simple( 1,     // number of dimensions of dataspace 
                                     dims,  // size of array in each dimension
                                     NULL); // allow stretching of data space (NULL=no)

    // period
    // create the dataset
    dataset_id = H5Dcreate( file_id,        // file identifier (or group identifier)
                            "/config/period",  // name of dataset (relative to group specified, if speficied)
                            H5T_NATIVE_LONG, // datatype to use when creating dataset
                            dataspace_id,   // dataspace identifier
                            H5P_DEFAULT,    // link creation property list (was dataset creating property list <=v1.6)
                            H5P_DEFAULT,    // dataset creation property list (added in HDF5v1.8)
                            H5P_DEFAULT);   // dataset access property list (added in HDF5v1.8)

    // write the dataset
    data2write_long[0]  = (long)period;
    status = H5Dwrite( dataset_id,          // dataset identifier
                       H5T_NATIVE_LONG,   // informs hdf about format of data in memory of computer
                       H5S_ALL,             // identifier of memory dataspace
                       H5S_ALL,             // file space identifier
                       H5P_DEFAULT,         // data transfer property list
                       data2write_long);     // pointer to data array
    // terminate access and free ressources/identifiers for dataset
    status = H5Dclose(dataset_id);

    // d_absorb
    dataset_id = H5Dcreate( file_id, "/config/d_absorb", H5T_NATIVE_LONG,
                            dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    data2write_long[0]  = (long)d_absorb;
    status = H5Dwrite( dataset_id, H5T_NATIVE_LONG,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data2write_long); 
    status = H5Dclose(dataset_id);

    // N_x
    dataset_id = H5Dcreate( file_id, "/config/N_x", H5T_NATIVE_LONG,
                            dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    data2write_long[0]  = (long)N_x;
    status = H5Dwrite( dataset_id, H5T_NATIVE_LONG,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data2write_long); 
    status = H5Dclose(dataset_id);

    // N_y
    dataset_id = H5Dcreate( file_id, "/config/N_y", H5T_NATIVE_LONG,
                            dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    data2write_long[0]  = (long)N_y;
    status = H5Dwrite( dataset_id, H5T_NATIVE_LONG,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data2write_long); 
    status = H5Dclose(dataset_id);

    // N_z
    dataset_id = H5Dcreate( file_id, "/config/N_z", H5T_NATIVE_LONG,
                            dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    data2write_long[0]  = (long)N_z;
    status = H5Dwrite( dataset_id, H5T_NATIVE_LONG,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data2write_long); 
    status = H5Dclose(dataset_id);


    // terminate access and free ressources/identifiers
    // data space
    status = H5Sclose(dataspace_id);
    // file 
    status = H5Fclose(file_id);
    
    return EXIT_SUCCESS;
}//#}}}
#endif


#ifdef HDF5
int readMyHDF( int dim0, int dim1, int dim2, char filename[], char dataset[], double array_3D[dim0][dim1][dim2]) {
    //#{{{

    // hdf handles
    hid_t           file_id, dset_id;
    herr_t          status;
    //hsize_t         dims[3] = { dim0, dim1, dim2};

    int             ii, jj;

    // open file using default properties
    file_id = H5Fopen( filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    // open dataset using default properties
    dset_id = H5Dopen( file_id, dataset, H5P_DEFAULT);

    // set the pointers to rows to the correct addresses
    // must be here, otherwise, data is not read correctly!
    // might be unnecessary, if matrix-allocation would be different (?)
    //for ( ii=1; ii<dims[0]; ii++ )
    //    array_2D[ii] = array_2D[0] + ii * dims[1];

    // Read the data using the default properties.
    status = H5Dread( dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      array_3D );

    // close the dataset
    status = H5Dclose( dset_id);

    // close the file
    status = H5Fclose( file_id);

    return EXIT_SUCCESS;
}//#}}}
#endif


int writeTimetraces2ascii( int dim0, int dim1, int t_end, double period, 
                           char filename[], double timetraces[dim0][dim1] ) {
//{{{

    size_t
        ii;

    FILE
        *file_pntr;

    // open file in w(rite) mode; might consider using a+ instead
    file_pntr   = fopen( filename, "w" );
    if (file_pntr == NULL) {
        printf( "ERROR: Unable to create file for timetraces.\n" );
        return EXIT_FAILURE;
    } else {
        // NOTE: if return value of printf < 0, then writing failed.
        //       might be good idea to implicetely check this
        //       e.g. if ( (fprintf( file_pntr, "a b c" )) < 0 ) ....
        fprintf( file_pntr, "# T  poynt_z1  poynt_z2  poynt_x1  poynt_x2  poynt_y1  poynt_y2  P_out\n" ); 
        for ( ii=0 ; ii<(t_end/(int)period) ; ++ii )
            fprintf( file_pntr, " %4d  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n",
                    (int)timetraces[ii][1], 
                    timetraces[ii][2], timetraces[ii][3],
                    timetraces[ii][4], timetraces[ii][5],
                    timetraces[ii][6], timetraces[ii][7],
                    (timetraces[ii][2]+timetraces[ii][3] + timetraces[ii][4]+timetraces[ii][5] + timetraces[ii][6]+timetraces[ii][7])
                  );
        if ((fclose(file_pntr)) == EOF) {
            printf( "ERROR: could not close file for timetraces.\n" );
        }
    }
    printf( "successfully written timetraces into %s\n", filename );
    return EXIT_SUCCESS;

}//}}}


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


