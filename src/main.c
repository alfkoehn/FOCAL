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
#include <time.h>


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


// setting boundary conditions, possible choices are
// 1: simple_abc
// 2: Mur
#define BOUNDARY 1

#define DETECTOR_ANTENNA_1D


#include "focal.h"
#include "antenna.h"
#include "grid_io.h"
#include "background_profiles.h"
#include "power_calc.h"
#include "initialize_grid.h"
#include "save_data.h"
#include "boundaries.h"
#include "antenna_detector.h"
#include "alloc-memory.h"


int main( int argc, char *argv[] ) {
//{{{

    clock_t start, end;
    double cpu_time_used;
    start = clock();

    /*Call structs*/
    struct gridConfiguration gridCfg;
    struct beamConfiguration beamCfg;
    struct namePath pathFile;
    struct antennaDetector ant_Detect;
    struct abcBoundary *ABC;

    int
        ii,jj,kk,
        t_int, T_wave, 

#ifdef _OPENMP
        n_threads,                          // number of threads that will be used (OpenMP)
#endif
        pwr_dect,

        len_str,                            // return value of snprintf
        opt_ret;                            // return value of getopt (reading input parameter)

    double

        poynt_x1, poynt_x2,
        poynt_y1, poynt_y2,
        poynt_z1, poynt_z2,
        poynt_z1_ref,

        power_abs_x1, power_abs_x2,
        power_abs_y1, power_abs_y2,
        power_abs_z1, power_abs_z2,
        power_abs_ref,

        omega_t;

    char
        dSet_name[PATH_MAX],
        filename_hdf5[PATH_MAX],            // filename of hdf5 file for output
        filename_timetraces[PATH_MAX],
        fullDir[2024];

    bool
        angle_zx_set,                       // is antAngle_zx set during call ?
        angle_zy_set;                       // is antAngle_zy set during call ?


    gridConfInit( &gridCfg, &pathFile, &beamCfg, &ant_Detect);  
    create_folder_path( &pathFile);
    initialize_antDetect( &ant_Detect, &gridCfg, &beamCfg);

    /*Allocate structs in memory*/
    if(pathFile.boundary == 1){
        ALLOC_1D(ABC, 1, abcBoundary);
    }
    else if(pathFile.boundary == 2){
        //ALLOC_1D(ABC, 1, abcBoundary);
    }
    
    setBoundary(&gridCfg, &pathFile, ABC);

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

    pwr_dect    = gridCfg.d_absorb;

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
    printf( "Boundary condition set to '%d'\n", pathFile.boundary );

    print_antDetect_info( &ant_Detect);

    //computeBoundary(  &gridCfg, PML, &pathFile);

#ifdef _OPENMP
#pragma omp parallel private(n_threads)
    {
    n_threads = omp_get_num_threads();
    printf( "number of threads that will be used (OpenMP) = %d\n", n_threads );
    }
#endif

    /*Initiate system evolution*/
    for (t_int=0 ; t_int <= gridCfg.t_end ; ++t_int) {
        
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
        if(pathFile.boundary == 1){
            apply_absorber(     &gridCfg, ABC->eco, EB_WAVE );
            apply_absorber_ref( &gridCfg, ABC->eco, EB_WAVE_ref );
        }

        // advance J
        advance_J( &gridCfg, EB_WAVE, J_B0, n_e );

        // advance B
        advance_B(     &gridCfg, EB_WAVE );
        advance_B_ref( &gridCfg, EB_WAVE_ref );
        
        // advance E
        advance_E(     &gridCfg, EB_WAVE,     J_B0 );
        advance_E_ref( &gridCfg, EB_WAVE_ref       );

        // optionally, apply numerical viscosity
        //apply_numerical_viscosity( &gridCfg, EB_WAVE );

        computeBoundary(  &gridCfg, &pathFile, ABC, EB_WAVE, t_int, EB_WAVE_ref);

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

        // store wavefields for detector antennas over the final 10 
        // oscillation periods, it was found previously that only one period
        // does not result in a too nice average
        /*if(ant_Detect.antDetect_1D == 1){
            
            if ( t_int >= (gridCfg.t_end-10*gridCfg.period) ) {
                if (ant_Detect.detAnt_01_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
                    detAnt1D_storeValues( &gridCfg, ant_Detect.detAnt_01_ypos, ant_Detect.detAnt_01_zpos,
                                        t_int,  
                                        EB_WAVE, detAnt_01_fields );
                }
                if (ant_Detect.detAnt_02_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
                    detAnt1D_storeValues( &gridCfg, ant_Detect.detAnt_01_ypos, ant_Detect.detAnt_02_zpos,
                                        t_int, 
                                        EB_WAVE, detAnt_02_fields );
                }
                if (ant_Detect.detAnt_03_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
                    detAnt1D_storeValues( &gridCfg, ant_Detect.detAnt_01_ypos, ant_Detect.detAnt_03_zpos,
                                        t_int,
                                        EB_WAVE, detAnt_03_fields );
                }
                if (ant_Detect.detAnt_04_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
                    detAnt1D_storeValues( &gridCfg, ant_Detect.detAnt_01_ypos, ant_Detect.detAnt_04_zpos,
                                        t_int,
                                        EB_WAVE, detAnt_04_fields );
                }
            }
        }*/
#ifdef DETECTOR_ANTENNA_1D
        // store wavefields for detector antennas over the final 10 
        // oscillation periods, it was found previously that only one period
        // does not result in a too nice average
        if ( t_int >= (gridCfg.t_end-10*gridCfg.period) ) {
            if (ant_Detect.detAnt_01_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
                detAnt1D_storeValues( &gridCfg, ant_Detect.detAnt_01_ypos, ant_Detect.detAnt_01_zpos,
                                      t_int,  
                                      EB_WAVE, detAnt_01_fields );
            }
            if (ant_Detect.detAnt_02_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
                detAnt1D_storeValues( &gridCfg, ant_Detect.detAnt_01_ypos, ant_Detect.detAnt_02_zpos,
                                      t_int, 
                                      EB_WAVE, detAnt_02_fields );
            }
            if (ant_Detect.detAnt_03_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
                detAnt1D_storeValues( &gridCfg, ant_Detect.detAnt_01_ypos, ant_Detect.detAnt_03_zpos,
                                      t_int,
                                      EB_WAVE, detAnt_03_fields );
            }
            if (ant_Detect.detAnt_04_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
                detAnt1D_storeValues( &gridCfg, ant_Detect.detAnt_01_ypos, ant_Detect.detAnt_04_zpos,
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

    if (snprintf(fullDir, sizeof(fullDir),"%s/%s", pathFile.projectPath, pathFile.foldername) >= sizeof(fullDir)) {
        fprintf(stderr, "Error: Directory path is too long.\n");
        return -1;
    }
    // write timetrace data into file
    // open file in w(rite) mode; might consider using a+ instead
    snprintf( filename_timetraces, sizeof(filename_timetraces), "%s/%s", fullDir, pathFile.file_trace);
    writeTimetraces2ascii( (gridCfg.t_end/(int)gridCfg.period), 8, gridCfg.t_end, gridCfg.period, 
                           filename_timetraces, timetraces );

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
    len_str = snprintf( filename_hdf5, sizeof(filename_hdf5), "%s/%s", fullDir, pathFile.file_hdf5);
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
    if (ant_Detect.detAnt_01_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
        detAnt1D_write2hdf5( gridCfg.Nx, filename_hdf5, "/detAnt_01" , 
                             ant_Detect.detAnt_01_ypos, ant_Detect.detAnt_01_zpos,
                             detAnt_01_fields );
    }
    if (ant_Detect.detAnt_02_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
        detAnt1D_write2hdf5( gridCfg.Nx, filename_hdf5, "/detAnt_02" , 
                             ant_Detect.detAnt_01_ypos, ant_Detect.detAnt_02_zpos,
                             detAnt_02_fields );
    }
    if (ant_Detect.detAnt_03_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
        detAnt1D_write2hdf5( gridCfg.Nx, filename_hdf5, "/detAnt_03" , 
                             ant_Detect.detAnt_01_ypos, ant_Detect.detAnt_03_zpos,
                             detAnt_03_fields );
    }
    if (ant_Detect.detAnt_04_zpos < (gridCfg.Nz - gridCfg.d_absorb)) {
        detAnt1D_write2hdf5( gridCfg.Nx, filename_hdf5, "/detAnt_04" , 
                             ant_Detect.detAnt_01_ypos, ant_Detect.detAnt_04_zpos,
                             detAnt_04_fields );
    }
#endif

    /*free structs from memory*/
    if(pathFile.boundary == 1){
        free( ABC );
    }
    else if(pathFile.boundary == 2){
        free( ABC );
    }

    free( EB_WAVE );
    printf( "freed EB_WAVE\n" );
    free( J_B0 );
    printf( "freed J_B0\n" );
    free( n_e );
    printf( "freed n_e\n" );
    free( data2save );
    printf( "freed data2save\n" );

    end = clock();
    cpu_time_used = ( (double)(end-start) )/CLOCKS_PER_SEC;
	printf("Program Running Time = %.2e s.\n", cpu_time_used);
    
    return EXIT_SUCCESS;
}//}}}


