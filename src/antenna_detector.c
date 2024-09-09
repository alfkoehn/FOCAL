#include "antenna_detector.h"

void initialize_antDetect(antennaDetector *ant_Detect, gridConfiguration *gridCfg, beamConfiguration *beamCfg){

    if(ant_Detect->antDetect_1D == 1){
    
        ant_Detect->detAnt_01_ypos  = beamCfg->ant_y;
        ant_Detect->detAnt_01_zpos  = beamCfg->ant_z+2;
        ant_Detect->detAnt_02_zpos  = round(beamCfg->ant_z+2 + 1*5*gridCfg->period); // steps of 5 cm for 28 GHz = 4.67*period
        ant_Detect->detAnt_03_zpos  = round(beamCfg->ant_z+2 + 2*5*gridCfg->period);
        ant_Detect->detAnt_04_zpos  = round(beamCfg->ant_z+2 + 3*5*gridCfg->period);
        // positions have to be even numbers, to ensure fields are accessed correctly
        if ((ant_Detect->detAnt_01_ypos % 2) != 0)  ++ant_Detect->detAnt_01_ypos;
        if ((ant_Detect->detAnt_01_zpos % 2) != 0)  ++ant_Detect->detAnt_01_zpos;
        if ((ant_Detect->detAnt_02_zpos % 2) != 0)  ++ant_Detect->detAnt_02_zpos;
        if ((ant_Detect->detAnt_03_zpos % 2) != 0)  ++ant_Detect->detAnt_03_zpos;
        if ((ant_Detect->detAnt_04_zpos % 2) != 0)  ++ant_Detect->detAnt_04_zpos;
        // issue a warning when detector antenna position is beyond Nz
        if (ant_Detect->detAnt_04_zpos > (gridCfg->Nz - gridCfg->d_absorb)) {
            printf( "ERROR: check the detector antenna positions into z direction\n" );
            printf( "       Nz-d_absorb = %d, detAnt_04_zpos = %d\n", 
                    gridCfg->Nz-gridCfg->d_absorb, ant_Detect->detAnt_04_zpos );

        }
    }

    if(ant_Detect->antDetect_1D != 1){
        printf( "Antenna Detector has not been initialized\n");
    }

}

void print_antDetect_info(antennaDetector *ant_Detect){

    if(ant_Detect->antDetect_1D == 1){
    printf( "detector antenna positions: z1 = %d, y1 = %d\n", ant_Detect->detAnt_01_zpos, ant_Detect->detAnt_01_ypos );
    printf( "detector antenna positions: z2 = %d, y1 = %d\n", ant_Detect->detAnt_02_zpos, ant_Detect->detAnt_01_ypos );
    printf( "detector antenna positions: z3 = %d, y1 = %d\n", ant_Detect->detAnt_03_zpos, ant_Detect->detAnt_01_ypos );
    printf( "detector antenna positions: z4 = %d, y1 = %d\n", ant_Detect->detAnt_04_zpos, ant_Detect->detAnt_01_ypos );
    }

}