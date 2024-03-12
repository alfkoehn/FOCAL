// header guard first to prevent multiple declarations
#ifndef FOCAL_H
#define FOCAL_H

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// define structures
typedef struct gridConfiguration {
    int
        Nx, Ny, Nz,     // maybe size_t would be better
        Nz_ref,
        d_absorb,
        t_end,
        ne_profile, B0_profile;
    double
        period,
        dx,dt;
} gridConfiguration;
typedef struct beamConfiguration {
    int
        exc_signal,
        ant_x, ant_y, ant_z,
        rampUpMethod;
    double
        antAngle_zy, antAngle_zx,
        ant_w0x, ant_w0y,
        z2waist;
} beamConfiguration;


int advance_J( gridConfiguration *gridCfg, 
               double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
               double J_B0[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz],
               double n_e[gridCfg->Nx/2][gridCfg->Ny/2][gridCfg->Nz/2] ); 

int advance_B( gridConfiguration *gridCfg, 
               double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] );

int advance_B_ref( gridConfiguration *gridCfg, 
                   double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref] );

int advance_E( gridConfiguration *gridCfg, 
               double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
               double J_B0[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] );

int advance_E_ref( gridConfiguration *gridCfg, 
                   double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref] ); 


#endif  // FOCAL_H

