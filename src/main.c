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
//PI values is defined in the math.h
#ifndef M_PI    
  #define M_PI 3.14159265358979323846   
#endif

//#define DETECTOR_ANTENNA_1D

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
#include "antenna_detector.h"

int main( int argc, char *argv[] ) {
//{{{

    gridConfiguration            *gridCfg;
    beamAntennaConfiguration     *beamCfg;
    saveData                     *saveDCfg;  
    antennaDetector              *antDetect; 
    boundaryVariables            *boundaryV;

    /*Alloc structs in memory*/
    ALLOC_1D( gridCfg, 1, gridConfiguration);
    ALLOC_1D( beamCfg, 1, beamAntennaConfiguration);
    ALLOC_1D( saveDCfg, 1, saveData);
    ALLOC_1D( antDetect, 1, antennaDetector );
    ALLOC_1D( boundaryV, 1, boundaryVariables);

    int
        t_int,  

#ifdef _OPENMP
        n_threads,                          // number of threads that will be used (OpenMP)
#endif
        pwr_dect;

    double

        poynt_x1, poynt_x2,
        poynt_y1, poynt_y2,
        poynt_z1, poynt_z2,
        poynt_z1_ref,

        power_abs_x1, power_abs_x2,
        power_abs_y1, power_abs_y2,
        power_abs_z1, power_abs_z2,
        power_abs_ref;


    // set-up grid
    control_init(  gridCfg, beamCfg, saveDCfg, antDetect );

    Y_at_X1     = .41;
    k0Ln_at_X1  = 6.;
    theta_at_X1 = 78.;

    create_folder( gridCfg, saveDCfg );
    init_boundary( gridCfg, boundaryV);
    init_antennaDetect( gridCfg, beamCfg, antDetect );

    // arrays realized as variable-length array (VLA)
    // E- and B-wavefield
    double (*EB_WAVE)[NY][NZ]           = calloc(NX, sizeof *EB_WAVE);
    double (*EB_WAVE_ref)[NY][NZ_REF]   = calloc(NX, sizeof *EB_WAVE_ref);
    // J-wavefield (in plasma) and background magnetic field
    double (*J_B0)[NY][NZ]              = calloc(NX, sizeof *J_B0);
    // background electron plasma density
    double (*n_e)[NY/2][NZ/2]           = calloc(NX/2, sizeof *n_e);
    // time traces
    double (*timetraces)[8]             = calloc((T_END/(int)period), sizeof *timetraces);

    pwr_dect    = d_absorb;

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

    init_antennaInjection( gridCfg, beamCfg );
    init_background_profiles(  gridCfg, beamCfg, n_e, J_B0 );

#ifdef _OPENMP
#pragma omp parallel private(n_threads)
    {
    n_threads = omp_get_num_threads();
    #pragma omp single
    printf( "number of threads that will be used (OpenMP) = %d\n", n_threads );
    }
#endif

    print_systemConfiguration( gridCfg, beamCfg );
    print_antennaDetec( antDetect );

    for (t_int=0 ; t_int <= T_END ; ++t_int) {
        
        control_antennaInjection(  gridCfg, beamCfg, t_int, EB_WAVE, EB_WAVE_ref ); //beam injection into grid
        advance_fields( gridCfg, EB_WAVE, EB_WAVE_ref, J_B0, n_e );                 //advance EM fields

        // optionally, apply numerical viscosity
        //apply_numerical_viscosity( &gridCfg, EB_WAVE );

        // apply absorbers
        advance_boundary(  gridCfg, boundaryV, EB_WAVE, EB_WAVE_ref );

        control_antennaDetect(  gridCfg, antDetect, t_int, EB_WAVE );

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

        save_field_toHDF5( gridCfg, saveDCfg, t_int, EB_WAVE ); //stores abs(E) into HDF5 file

    } // end of time loop
  
    // write timetrace data into file
    writeConsole_timetraces( (T_END/(int)period), col_for_timetraces, T_END, period, timetraces );
    control_save( gridCfg, beamCfg ,saveDCfg, timetraces, n_e, J_B0 );

    save_AntDetect( gridCfg, saveDCfg, antDetect );

    free_boundary( gridCfg );
    free_antDetect( gridCfg, antDetect );
    free( EB_WAVE );
    printf( "freed EB_WAVE\n" );
    free( J_B0 );
    printf( "freed J_B0\n" );
    free( n_e );
    printf( "freed n_e\n" );
    
    return EXIT_SUCCESS;
}//}}}
