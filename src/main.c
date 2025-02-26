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

//Include C libraries
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

//Include FOCAL modules headers
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

    //Call structures
    gridConfiguration            *gridCfg;
    beamAntennaConfiguration     *beamCfg;
    saveData                     *saveDCfg;  
    antennaDetector              *antDetect; 
    boundaryVariables            *boundaryV;
    powerValues                  *powerVal;

    /*Alloc structs in memory*/
    ALLOC_1D( gridCfg, 1, gridConfiguration);
    ALLOC_1D( beamCfg, 1, beamAntennaConfiguration);
    ALLOC_1D( saveDCfg, 1, saveData);
    ALLOC_1D( antDetect, 1, antennaDetector );
    ALLOC_1D( boundaryV, 1, boundaryVariables);
    ALLOC_1D( powerVal, 1, powerValues);

    int  
#ifdef _OPENMP
        n_threads,                          // number of threads that will be used (OpenMP)
#endif
        t_int;

    // set-up grid (read values from JSON)
    control_init(  gridCfg, beamCfg, saveDCfg, antDetect );     //function in INIT_MODULE.C
    create_folder( gridCfg, saveDCfg );                         //function in SAVE_DATA.C

    printf("----------------Initializing Profiles---------------\n");
    init_boundary( gridCfg, boundaryV);                         //function in BOUNDARY_MODULE.C
    init_antennaDetect( gridCfg, beamCfg, antDetect );          //function in ANTENNA_DETECTOR.C
    init_powerValues( gridCfg, powerVal );                      //function in POWER_CALC.C      
    init_antennaInjection( gridCfg, beamCfg );                  //function in ANTENNA.C

    // arrays realized as variable-length array (VLA)
    // E- and B-wavefield
    double (*EB_WAVE)[NY][NZ]           = calloc(NX, sizeof *EB_WAVE);
    double (*EB_WAVE_ref)[NY][NZ_REF]   = calloc(NX, sizeof *EB_WAVE_ref);
    // J-wavefield (in plasma) and background magnetic field
    double (*J_B0)[NY][NZ]              = calloc(NX, sizeof *J_B0);
    // background electron plasma density
    double (*n_e)[NY/2][NZ/2]           = calloc(NX/2, sizeof *n_e);

    init_background_profiles(  gridCfg, beamCfg, saveDCfg, n_e, J_B0 );   //function in BACKGROUND_PROFILES.C

    //Simulation values print to terminal
    print_systemConfiguration( gridCfg, beamCfg, powerVal );    //function in INIT_MODULE.C
    print_antennaDetec( antDetect );                            //function in ANTENNA_DETECTOR.C

    printf("--------------------Simulation--------------------\n");

#ifdef _OPENMP
#pragma omp parallel private(n_threads)
    {
    n_threads = omp_get_num_threads();
    #pragma omp single
    printf( "number of threads that will be used (OpenMP) = %d\n", n_threads );
    }
#endif

    //System's time evolution
    for ( t_int=0 ; t_int <= T_END ; ++t_int ) {
        
        //Beam injection to grid
        control_antennaInjection(  gridCfg, beamCfg, t_int, EB_WAVE, EB_WAVE_ref ); //function in ANTENNA.C
        advance_fields( gridCfg, EB_WAVE, EB_WAVE_ref, J_B0, n_e );                 //advance EM fields. function in FOCAL.C

        //optionally, apply numerical viscosity
        //apply_numerical_viscosity( &gridCfg, EB_WAVE );

        advance_boundary(  gridCfg, boundaryV, EB_WAVE, EB_WAVE_ref );              //function in BOUNDARY_MODULE.C
        control_antennaDetect(  gridCfg, antDetect, t_int, EB_WAVE );               //function in ANTENNA_DETECTOR.C
        compute_power( gridCfg, beamCfg, powerVal, t_int, EB_WAVE, EB_WAVE_ref );   //function in POWER_CALC.C
        //stores abs(E) into HDF5 file
        save_field_toHDF5( gridCfg, saveDCfg, t_int, EB_WAVE );                     //function in SAVE_DATA.C

    } // end of time loop
    
    printf("--------------------Finished!-------------------\n");

    // write timetrace data into file
    write_timetraces( gridCfg, saveDCfg );                                  //function in POWER_CALC.C
    save_SimData( gridCfg, beamCfg ,saveDCfg, n_e, J_B0 );                  //function in SAVE_DATA.C      
    save_AntDetect( gridCfg, saveDCfg, antDetect );                         //function in ANTENNA_DETECTOR.C

    //free allocated arrays
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
