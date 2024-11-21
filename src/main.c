/**
 * Author:      Alf Köhn-Seemann
 * Email:       koehn@igvp.uni-stuttgart.de
 * Copyright:   University of Stuttgart
 * 
 * This is a 3D FDTD code for simulating electromagnetic waves in cold 
 * magnetized plasmas. The code is called FOCAL.
 *
 * NOTE: FOCAL is still in development (but physicswise delivering 
 *       correct, i.e. benchmarked in some cases, results).
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


// setting boundary conditions, possible choices are
// 1: simple_abc
// 2: Mur
//#define BOUNDARY 2

#define DETECTOR_ANTENNA_1D

#include "focal-struct.h"
#include "macros-grid.h"
#include "alloc-memory.h"
#include "init_module.h"
#include "focal.h"
#include "antenna.h"
#include "grid_io.h"
#include "background_profiles.h"
#include "power_calc.h"
#include "save_data.h"
#include "boundary_module.h"

int main( int argc, char *argv[] ) {
//{{{

    gridConfiguration            *gridCfg;
    beamAntennaConfiguration     *beamCfg;
    saveData                     *saveDCfg;  
    antennaDetector              *antDetect; 
    boundaryVariables            *boundaryV;
    /*struct codeDiagnostics              *diagnostic;*/

    /*Alloc structs in memory*/
    ALLOC_1D( gridCfg, 1, gridConfiguration);
    ALLOC_1D( beamCfg, 1, beamAntennaConfiguration);
    ALLOC_1D( saveDCfg, 1, saveData);
    ALLOC_1D( antDetect, 1, antennaDetector );
    ALLOC_1D( boundaryV, 1, boundaryVariables);
    /*ALLOC_1D( diagnostic, 1, codeDiagnostics );*/

    int
        t_int, T_wave, 

#ifdef _OPENMP
        n_threads,                          // number of threads that will be used (OpenMP)
#endif
        pwr_dect,

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

    bool
        angle_zx_set,                       // is antAngle_zx set during call ?
        angle_zy_set;                       // is antAngle_zy set during call ?

    // set-up grid
    control_init(  gridCfg, beamCfg, saveDCfg, antDetect );

    Y_at_X1     = .41;
    k0Ln_at_X1  = 6.;
    theta_at_X1 = 78.;

    create_folder( gridCfg, saveDCfg );
    init_boundary( gridCfg, boundaryV);

    // arrays realized as variable-length array (VLA)
    // E- and B-wavefield
    double (*EB_WAVE)[NY][NZ]           = calloc(NX, sizeof *EB_WAVE);
    double (*EB_WAVE_ref)[NY][NZ_REF]   = calloc(NX, sizeof *EB_WAVE_ref);
    // J-wavefield (in plasma) and background magnetic field
    double (*J_B0)[NY][NZ]              = calloc(NX, sizeof *J_B0);
    // background electron plasma density
    double (*n_e)[NY/2][NZ/2]           = calloc(NX/2, sizeof *n_e);
    // used when writing data into hdf5-files
    //double (*data2save)[NY/2][NZ/2]     = calloc(NX/2, sizeof *data2save);
    // antenna: envelope of injected field
    double (*antField_xy)[NY/2]         = calloc(NX/2, sizeof *antField_xy);
    // antenna: phase terms 
    double (*antPhaseTerms)[NY/2]       = calloc(NX/2, sizeof *antPhaseTerms);
    // time traces
    double (*timetraces)[8]             = calloc((T_END/(int)period), sizeof *timetraces);

    // old E-fields required for Mur's boundary condition
/*#if BOUNDARY == 2
    double (*E_Xdir_OLD)[NY][NZ]            = calloc(8,  sizeof *E_Xdir_OLD);
    double (*E_Ydir_OLD)[8][NZ]             = calloc(NX, sizeof *E_Ydir_OLD);
    double (*E_Zdir_OLD)[NY][8]             = calloc(NX, sizeof *E_Zdir_OLD);
    double (*E_Xdir_OLD_ref)[NY][NZ_REF]    = calloc(8,  sizeof *E_Xdir_OLD_ref);
    double (*E_Ydir_OLD_ref)[8][NZ_REF]     = calloc(NX, sizeof *E_Ydir_OLD_ref);
    double (*E_Zdir_OLD_ref)[NY][8]         = calloc(NX, sizeof *E_Zdir_OLD_ref);
#endif*/

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
            case 'a': antAngle_zx   = atof(optarg);
                      angle_zx_set  = true;
                      break;
            case 'b': antAngle_zy   = atof(optarg);
                      angle_zy_set  = true;
                      break;
        }
    }
    if ( argc > 1 ) {
        printf( "following parameters were set during call: \n" );
        if (angle_zx_set)   printf( "    antAngle_zx = %f\n", antAngle_zx );
        if (angle_zy_set)   printf( "    antAngle_zy = %f\n", antAngle_zy );
    }

    pwr_dect    = d_absorb;

#ifdef DETECTOR_ANTENNA_1D
    detAnt_01_ypos  = ant_y;
    detAnt_01_zpos  = ant_z+2;
    detAnt_02_zpos  = round(ant_z+2 + 1*5*period); // steps of 5 cm for 28 GHz = 4.67*period
    detAnt_03_zpos  = round(ant_z+2 + 2*5*period);
    detAnt_04_zpos  = round(ant_z+2 + 3*5*period);
    // positions have to be even numbers, to ensure fields are accessed correctly
    if ((detAnt_01_ypos % 2) != 0)  ++detAnt_01_ypos;
    if ((detAnt_01_zpos % 2) != 0)  ++detAnt_01_zpos;
    if ((detAnt_02_zpos % 2) != 0)  ++detAnt_02_zpos;
    if ((detAnt_03_zpos % 2) != 0)  ++detAnt_03_zpos;
    if ((detAnt_04_zpos % 2) != 0)  ++detAnt_04_zpos;
    // issue a warning when detector antenna position is beyond NZ
    if (detAnt_04_zpos > (NZ - d_absorb)) {
        printf( "ERROR: check the detector antenna positions into z direction\n" );
        printf( "       NZ-d_absorb = %d, detAnt_04_zpos = %d\n", 
                NZ - d_absorb, detAnt_04_zpos );
    }
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
    printf( "...done setting all variables to 0\n" );

   printf( "starting do define antenna field...\n" );
    make_antenna_profile( gridCfg, beamCfg, 
                          antField_xy, antPhaseTerms );
    printf( "...done defining antenna field\n" );

    printf( "starting defining background plasma density\n" );
            // ne_profile: 1 = plasma mirror
            //             2 = linearly increasing profile
    make_density_profile( gridCfg,  
            // cntrl_para: ne_profile=1 --> 0: plane mirror; oblique mirror: -.36397; 20 degrees: -.17633
            //             ne_profile=2 --> k0*Ln: 25
            k0Ln_at_X1,
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
            gridCfg,
            // cntrl_para: B0_profile=1 --> value of Y
            Y_at_X1, 
            J_B0 );
    printf( "...done defining background magnetic field\n" );

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

    print_systemConfiguration( gridCfg, beamCfg );

    for (t_int=0 ; t_int <= T_END ; ++t_int) {
        
        omega_t += 2.*M_PI/period;

        // to avoid precision problems when a lot of pi's are summed up        
        if (omega_t >= 2.*M_PI) {
            omega_t    += -2.*M_PI;
            T_wave     += 1;
            //printf("status: number of oscillation periods: %d (t_int= %d) \n",T_wave,t_int);
        }

        // add source
        add_source( gridCfg, beamCfg,
                    t_int, omega_t, 
                    antField_xy, antPhaseTerms, EB_WAVE );
        add_source_ref( gridCfg, beamCfg,
                        t_int, omega_t, 
                        antField_xy, antPhaseTerms, EB_WAVE_ref );


        advance_fields( gridCfg, EB_WAVE, EB_WAVE_ref, J_B0, n_e );

        // optionally, apply numerical viscosity
        //apply_numerical_viscosity( &gridCfg, EB_WAVE );

        // apply Mur's boundary conditions
/*#if BOUNDARY == 2
        abc_Mur_1st( gridCfg, "x1x2y1y2z1z2",  
                     EB_WAVE, E_Xdir_OLD, E_Ydir_OLD, E_Zdir_OLD );
        abc_Mur_1st_ref( gridCfg, 
                         EB_WAVE_ref, E_Xdir_OLD_ref, E_Ydir_OLD_ref, E_Zdir_OLD_ref );
        abc_Mur_saveOldE_xdir(    gridCfg, EB_WAVE, E_Xdir_OLD );
        abc_Mur_saveOldE_ydir(    gridCfg, EB_WAVE, E_Ydir_OLD );
        abc_Mur_saveOldE_zdir(    gridCfg, EB_WAVE, E_Zdir_OLD );
        abc_Mur_saveOldEref_xdir( gridCfg, EB_WAVE_ref, E_Xdir_OLD_ref );
        abc_Mur_saveOldEref_ydir( gridCfg, EB_WAVE_ref, E_Ydir_OLD_ref );
        abc_Mur_saveOldEref_zdir( gridCfg, EB_WAVE_ref, E_Zdir_OLD_ref );
#endif*/

        // apply absorbers
        advance_boundary(  gridCfg, boundaryV, EB_WAVE, EB_WAVE_ref );

#ifdef DETECTOR_ANTENNA_1D
        // store wavefields for detector antennas over the final 10 
        // oscillation periods, it was found previously that only one period
        // does not result in a too nice average
        if ( t_int >= ( T_END - 10*period ) ) {
            if (detAnt_01_zpos < ( NZ - d_absorb)) {
                detAnt1D_storeValues( gridCfg, detAnt_01_ypos, detAnt_01_zpos,
                                      t_int,  
                                      EB_WAVE, detAnt_01_fields );
            }
            if (detAnt_02_zpos < ( NZ - d_absorb)) {
                detAnt1D_storeValues( gridCfg, detAnt_01_ypos, detAnt_02_zpos,
                                      t_int, 
                                      EB_WAVE, detAnt_02_fields );
            }
            if (detAnt_03_zpos < ( NZ - d_absorb)) {
                detAnt1D_storeValues( gridCfg, detAnt_01_ypos, detAnt_03_zpos,
                                      t_int,
                                      EB_WAVE, detAnt_03_fields );
            }
            if (detAnt_04_zpos < ( NZ - d_absorb)) {
                detAnt1D_storeValues( gridCfg, detAnt_01_ypos, detAnt_04_zpos,
                                      t_int,
                                      EB_WAVE, detAnt_04_fields );
            }
        }
#endif

        // IQ detector for power detection
        if ( t_int >= 20*period ) {
            // z1-plane and z2-plane
            poynt_z1_ref    = calc_poynt_4( gridCfg, pwr_dect, "ref_z1", EB_WAVE, EB_WAVE_ref );
            poynt_z1        = calc_poynt_4( gridCfg, pwr_dect, "z1",     EB_WAVE, EB_WAVE_ref );
            poynt_z2        = calc_poynt_4( gridCfg, pwr_dect, "z2",     EB_WAVE, EB_WAVE_ref );
            // x1-plane and x2-plane
            poynt_x1        = calc_poynt_4( gridCfg, pwr_dect, "x1", EB_WAVE, EB_WAVE_ref );
            poynt_x2        = calc_poynt_4( gridCfg, pwr_dect, "x2", EB_WAVE, EB_WAVE_ref );
            // y1-plane and y2-plane
            poynt_y1        = calc_poynt_4( gridCfg, pwr_dect, "y1", EB_WAVE, EB_WAVE_ref );
            poynt_y2        = calc_poynt_4( gridCfg, pwr_dect, "y2", EB_WAVE, EB_WAVE_ref );

            
            //printf( "t = %d, power_abs_ref = %13.5e, power_abs_z1 = %13.5e, power_abs_z2 = %13.5e, poynt_z1 = %13.5e, poynt_z2 = %13.5e\n",
            //        t_int, power_abs_ref, power_abs_z1, power_abs_z2, poynt_z1, poynt_z2 );

            power_abs_ref   = .99*power_abs_ref + .01*poynt_z1_ref;
            power_abs_z1    = .99*power_abs_z1  + .01*poynt_z1;
            power_abs_z2    = .99*power_abs_z2  + .01*poynt_z2;
            power_abs_x1    = .99*power_abs_x1  + .01*poynt_x1;
            power_abs_x2    = .99*power_abs_x2  + .01*poynt_x2;
            power_abs_y1    = .99*power_abs_y1  + .01*poynt_y1;
            power_abs_y2    = .99*power_abs_y2  + .01*poynt_y2;

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
            timetraces[T_wave][0]   = (double)t_int;
            timetraces[T_wave][1]   = (double)T_wave;
            timetraces[T_wave][2]   = power_abs_z1/power_abs_ref;
            timetraces[T_wave][3]   = power_abs_z2/power_abs_ref;
            timetraces[T_wave][4]   = power_abs_x1/power_abs_ref;
            timetraces[T_wave][5]   = power_abs_x2/power_abs_ref;
            timetraces[T_wave][6]   = power_abs_y1/power_abs_ref;
            timetraces[T_wave][7]   = power_abs_y2/power_abs_ref;

        }

        save_field_toHDF5( gridCfg, saveDCfg, t_int, EB_WAVE );
        //writeUPMLdata( gridCfg, saveDCfg, EB_WAVE, t_int );

    } // end of time loop

    
    // write timetrace data into file
    writeConsole_timetraces( (T_END/(int)period), col_for_timetraces, T_END, period, timetraces );
    control_save( gridCfg, beamCfg ,saveDCfg, antDetect, timetraces, n_e, J_B0,
                  detAnt_01_fields, detAnt_02_fields, detAnt_03_fields, detAnt_04_fields );


    free( EB_WAVE );
    printf( "freed EB_WAVE\n" );
    free( J_B0 );
    printf( "freed J_B0\n" );
    free( n_e );
    printf( "freed n_e\n" );

    free_boundary( gridCfg );
    
    return EXIT_SUCCESS;
}//}}}


