#ifndef MACROS_BOUNDARY_H
#define MACROS_BOUNDARY_H

#include "focal-struct.h"

/*Macros for grid configuration*/
#define NxG(gridCfg)            gridCfg->Nx
#define NyG(gridCfg)            gridCfg->Ny
#define NzG(gridCfg)            gridCfg->Nz
#define Nz_refG(gridCfg)        gridCfg->Nz_ref
#define d_absorbG(gridCfg)      gridCfg->d_absorb
#define t_endG(gridCfg)         gridCfg->t_end
#define periodG(gridCfg)        gridCfg->period
#define dxG(gridCfg)            gridCfg->dx
#define dtG(gridCfg)            gridCfg->dt
#define ne_profileG(gridCfg)    gridCfg->ne_profile
#define ne_0G(gridCfg)          gridCfg->ne_0
#define B0_profileG(gridCfg)    gridCfg->B0_profile
#define B0_valueG(gridCfg)      gridCfg->B0_value
#define boundaryG(gridCfg)      gridCfg->sel_boundary

#define NX                      NxG(gridCfg)            
#define NY                      NyG(gridCfg)          
#define NZ                      NzG(gridCfg)
#define NZ_REF                  Nz_refG(gridCfg)  
#define d_absorb                d_absorbG(gridCfg)
#define T_END                   t_endG(gridCfg)
#define period                  periodG(gridCfg)
#define DX                      dxG(gridCfg)
#define DT                      dtG(gridCfg)
#define ne_profile              ne_profileG(gridCfg)
#define ne_0                    ne_0G(gridCfg)
#define B0_profile              B0_profileG(gridCfg)
#define B0_value                B0_valueG(gridCfg)
#define sel_boundary            boundaryG(gridCfg)

/*Macros for save data*/
#define projectPathSt(saveDCfg)         saveDCfg->projectPath
#define foldernameSt(saveDCfg)          saveDCfg->foldername
#define file_hdf5St(saveDCfg)           saveDCfg->file_hdf5
#define file_traceSt(saveDCfg)          saveDCfg->file_trace
#define t_saveSt(saveDCfg)              saveDCfg->t_save
#define col_for_timetracesSt(saveDCfg)  saveDCfg->col_for_timetraces

#define projectPath                     projectPathSt(saveDCfg)
#define foldername                      foldernameSt(saveDCfg)
#define file_hdf5                       file_hdf5St(saveDCfg)
#define file_trace                      file_traceSt(saveDCfg)
#define t_save                          t_saveSt(saveDCfg)
#define col_for_timetraces              col_for_timetracesSt(saveDCfg)

/*Macros for antenna injection*/
/*#define antFieldBG(beamCfg,i,j)         beamCfg->antField_xy[( (i) * (NY/2) ) + j ]
#define antPhaseBG(beamCfg,i,j)         beamCfg->antPhaseTerms[( (i) * (NY/2) ) + j ]
#define T_waveBG(beamCfg)               beamCfg->T_wave*/
#define exc_signalBG(beamCfg)           beamCfg->exc_signal
#define ant_xBG(beamCfg)                beamCfg->ant_x
#define ant_yBG(beamCfg)                beamCfg->ant_y
#define ant_zBG(beamCfg)                beamCfg->ant_z
#define rampUpMBG(beamCfg)              beamCfg->rampUpMethod
//#define omega_tBG(beamCfg)              beamCfg->omega_t
#define antAngle_zxBG(beamCfg)          beamCfg->antAngle_zx
#define antAngle_zyBG(beamCfg)          beamCfg->antAngle_zy
#define ant_w0xBG(beamCfg)              beamCfg->ant_w0x
#define ant_w0yBG(beamCfg)              beamCfg->ant_w0y
#define z2waistBG(beamCfg)              beamCfg->z2waist
#define Y_at_X1BG(beamCfg)              beamCfg->Y_at_X1
#define k0Ln_at_X1BG(beamCfg)           beamCfg->k0Ln_at_X1
#define theta_at_X1BG(beamCfg)          beamCfg->theta_at_X1

/*#define antField_xy(i,j)                antFieldBG(beamCfg,i,j)         
#define antPhaseTerms(i,j)              antPhaseBG(beamCfg,i,j)         
#define T_wave                          T_waveBG(beamCfg)*/               
#define exc_signal                      exc_signalBG(beamCfg)           
#define ant_x                           ant_xBG(beamCfg)               
#define ant_y                           ant_yBG(beamCfg)               
#define ant_z                           ant_zBG(beamCfg)               
#define rampUpMethod                    rampUpMBG(beamCfg)              
//#define t_omega                         omega_tBG(beamCfg)              
#define antAngle_zx                     antAngle_zxBG(beamCfg)          
#define antAngle_zy                     antAngle_zyBG(beamCfg)          
#define ant_w0x                         ant_w0xBG(beamCfg)              
#define ant_w0y                         ant_w0yBG(beamCfg)              
#define z2waist                         z2waistBG(beamCfg)  
#define Y_at_X1                         Y_at_X1BG(beamCfg)              
#define k0Ln_at_X1                      k0Ln_at_X1BG(beamCfg)           
#define theta_at_X1                     theta_at_X1BG(beamCfg)                     

/*Macros for ABC*/
#define ecoBV(boundaryV)                        boundaryV->eco

#define eco                                     ecoBV(boundaryV)

/*Macros for Mur boundary*/
/*#define E_Xdir_OLD_G(boundaryG,i,j,k)           boundaryG->E_Xdir_OLD[((i) * (NY) + j) * (NZ) + k]
#define E_Ydir_OLD_G(boundaryG,i,j,k)           boundaryG->E_Ydir_OLD[((i) * (d_absorb) + j) * (NZ) + k]
#define E_Zdir_OLD_G(boundaryG,i,j,k)           boundaryG->E_Zdir_OLD[((i) * (NY) + j) * (d_absorb) + k]
#define E_Xdir_OLD_ref_G(boundaryG,i,j,k)       boundaryG->E_Xdir_OLD_ref[((i) * (NY) + j) * (NZ_REF) + k]
#define E_Ydir_OLD_ref_G(boundaryG,i,j,k)       boundaryG->E_Ydir_OLD_ref[((i) * (d_absorb) + j) * (NZ_REF) + k]
#define E_Zdir_OLD_ref_G(boundaryG,i,j,k)       boundaryG->E_Zdir_OLD_ref[((i) * (NY) + j) * (d_absorb) + k]

#define E_Xdir_OLD(i,j,k)                       E_Xdir_OLD_G(boundaryG,i,j,k)
#define E_Ydir_OLD(i,j,k)                       E_Ydir_OLD_G(boundaryG,i,j,k)
#define E_Zdir_OLD(i,j,k)                       E_Zdir_OLD_G(boundaryG,i,j,k)
#define E_Xdir_OLD_ref(i,j,k)                   E_Xdir_OLD_ref_G(boundaryG,i,j,k)
#define E_Ydir_OLD_ref(i,j,k)                   E_Ydir_OLD_ref_G(boundaryG,i,j,k)
#define E_Zdir_OLD_ref(i,j,k)                   E_Zdir_OLD_ref_G(boundaryG,i,j,k)*/

/*Macros for UPML boundary layer*/
/*#define DH_WAVEstr(PMLG,i,j,k)                  boundaryG->DH_WAVE[((i) * (NY) + j) * (NZ) + k]
#define DH_WAVE_refStr(PMLG,i,j,k)              boundaryG->DH_WAVE_ref[((i) * (NY) + j) * (NZ_REF) + k]*/
#define F1xStr(boundaryV,i)                     boundaryV->F1x[i]
#define F1yStr(boundaryV,j)                     boundaryV->F1y[j]
#define F1zStr(boundaryV,k)                     boundaryV->F1z[k]
#define F2xStr(boundaryV,i)                     boundaryV->F2x[i]
#define F2yStr(boundaryV,j)                     boundaryV->F2y[j]
#define F2zStr(boundaryV,k)                     boundaryV->F2z[k]
#define CxStr(boundaryV,i)                      boundaryV->Cx[i]
#define CyStr(boundaryV,j)                      boundaryV->Cy[j]
#define CzStr(boundaryV,k)                      boundaryV->Cz[k]

#define F1zrStr(boundaryV,k)                    boundaryV->F1zr[k]
#define F2zrStr(boundaryV,k)                    boundaryV->F2zr[k]
#define CzrStr(boundaryV,k)                     boundaryV->Czr[k]

/*#define DH_WAVE(i,j,k)                          DH_WAVEstr(PMLG,i,j,k)
#define DH_WAVE_ref(i,j,k)                      DH_WAVE_refStr(PMLG,i,j,k)  */            
#define F1x(i)                                  F1xStr(boundaryV,i) 
#define F1y(j)                                  F1yStr(boundaryV,j)
#define F1z(k)                                  F1zStr(boundaryV,k)
#define F2x(i)                                  F2xStr(boundaryV,i)
#define F2y(j)                                  F2yStr(boundaryV,j) 
#define F2z(k)                                  F2zStr(boundaryV,k)
#define Cx(i)                                   CxStr(boundaryV,i)
#define Cy(j)                                   CyStr(boundaryV,j)
#define Cz(k)                                   CzStr(boundaryV,k)

#define F1zr(k)                                 F1zrStr(boundaryV,k)
#define F2zr(k)                                 F2zrStr(boundaryV,k)
#define Czr(k)                                  CzrStr(boundaryV,k)


/*Macros for Antenna detector*/
/*#define DET_ANT_ACCES(antDetect, id, i, j)  \
    ((id == FIELD_01) ? ( detAnt_01_fields ) : \
     (id == FIELD_02) ? ( detAnt_01_fields ) : \
     (id == FIELD_03) ? ( detAnt_03_fields ) : \
     (id == FIELD_04) ? ( detAnt_03_fields ) : 0 )*/

/*#define antDetect01EBG(antDetect,i,j)           antDetect->detAnt_01_fields[ ((i) * 5) + j ]
#define antDetect02EBG(antDetect,i,j)           antDetect->detAnt_02_fields[ ((i) * 5) + j ]
#define antDetect03EBG(antDetect,i,j)           antDetect->detAnt_03_fields[ ((i) * 5) + j ]
#define antDetect04EBG(antDetect,i,j)           antDetect->detAnt_04_fields[ ((i) * 5) + j ]*/
#define antDetect_1DG(antDetect)                antDetect->antDetect_1D
#define detAnt01zG(antDetect)                   antDetect->detAnt_01_zpos
#define detAnt02zG(antDetect)                   antDetect->detAnt_02_zpos
#define detAnt03zG(antDetect)                   antDetect->detAnt_03_zpos
#define detAnt04zG(antDetect)                   antDetect->detAnt_04_zpos
#define detAnt01yG(antDetect)                   antDetect->detAnt_01_ypos

/*#define detAnt_01_Fields(i,j)                   antDetect01EBG(antDetect,i,j)
#define detAnt_02_Fields(i,j)                   antDetect02EBG(antDetect,i,j)
#define detAnt_03_Fields(i,j)                   antDetect03EBG(antDetect,i,j)           
#define detAnt_04_Fields(i,j)                   antDetect04EBG(antDetect,i,j)*/
#define antDetect_1D                            antDetect_1DG(antDetect)
#define detAnt_01_zpos                          detAnt01zG(antDetect)
#define detAnt_02_zpos                          detAnt02zG(antDetect)
#define detAnt_03_zpos                          detAnt03zG(antDetect)
#define detAnt_04_zpos                          detAnt04zG(antDetect)
#define detAnt_01_ypos                          detAnt01yG(antDetect)


/*Macros for power struct*/
/*#define powerDectS(powerValStr)                 powerValStr->pwr_dect
#define powerAbsX1S(powerValStr)                powerValStr->power_abs_x1
#define powerAbsX2S(powerValStr)                powerValStr->power_abs_x2
#define powerAbsY1S(powerValStr)                powerValStr->power_abs_y1
#define powerAbsY2S(powerValStr)                powerValStr->power_abs_y2
#define powerAbsZ1S(powerValStr)                powerValStr->power_abs_z1
#define powerAbsZ2S(powerValStr)                powerValStr->power_abs_z2
#define powerAbsRefS(powerValStr)               powerValStr->power_abs_ref
#define powerPoyX1S(powerValStr)                powerValStr->poynt_x1
#define powerPoyX2S(powerValStr)                powerValStr->poynt_x2
#define powerPoyY1S(powerValStr)                powerValStr->poynt_y1
#define powerPoyY2S(powerValStr)                powerValStr->poynt_y2
#define powerPoyZ1S(powerValStr)                powerValStr->poynt_z1
#define powerPoyZRS(powerValStr)                powerValStr->poynt_z1_ref
#define powerPoyZ2S(powerValStr)                powerValStr->poynt_z2

#define pwr_dect                                powerDectS(powerValStr)
#define power_abs_x1                            powerAbsX1S(powerValStr)
#define power_abs_x2                            powerAbsX2S(powerValStr)                
#define power_abs_y1                            powerAbsY1S(powerValStr)               
#define power_abs_y2                            powerAbsY2S(powerValStr)                
#define power_abs_z1                            powerAbsZ1S(powerValStr)                
#define power_abs_z2                            powerAbsZ2S(powerValStr)                
#define power_abs_ref                           powerAbsRefS(powerValStr)               
#define poynt_x1                                powerPoyX1S(powerValStr)                
#define poynt_x2                                powerPoyX2S(powerValStr)                
#define poynt_y1                                powerPoyY1S(powerValStr)
#define poynt_y2                                powerPoyY2S(powerValStr)                
#define poynt_z1                                powerPoyZ1S(powerValStr)   
#define poynt_z1_ref                            powerPoyZRS(powerValStr) 
#define poynt_z2                                powerPoyZ2S(powerValStr) */ 

/*Diagnostics Structure*/
/*#define TotalEstr(diagnostic, i)                diagnostic->E[i]

#define E(i)                                    TotalEstr(diagnostic, i)*/

#endif
