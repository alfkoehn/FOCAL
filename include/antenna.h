#ifndef ANTENNA_H
#define ANTANNA_H

int make_antenna_profile( gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg,
                          double antField_xy[gridCfg->Nx/2][gridCfg->Ny/2], double antPhaseTerms[gridCfg->Nx/2][gridCfg->Ny/2] );

int add_source( gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg, 
                int t_int, double omega_t, 
                double antField_xy[gridCfg->Nx/2][gridCfg->Ny/2], 
                double antPhaseTerms[gridCfg->Nx/2][gridCfg->Ny/2],
                double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] );

int add_source_ref( gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg, 
                    int t_int, double omega_t, 
                    double antField_xy[gridCfg->Nx/2][gridCfg->Ny/2], 
                    double antPhaseTerms[gridCfg->Nx/2][gridCfg->Ny/2],
                    double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref] );

double antenna_field_rampup( int rampUpMethod, double period, int t_int );

double antenna_calcHansenExEy_O( double theta_rad, double Y );

#endif  // ANTENNA_H
