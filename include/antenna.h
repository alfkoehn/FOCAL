#ifndef ANTENNA_H
#define ANTANNA_H

int make_antenna_profile( gridConfiguration *gridCfg, beamConfiguration *beamCfg,
                          double antField_xy[gridCfg->Nx/2][gridCfg->Ny/2], double antPhaseTerms[gridCfg->Nx/2][gridCfg->Ny/2] );

#endif  // ANTENNA_H
