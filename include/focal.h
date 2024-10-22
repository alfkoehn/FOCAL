// header guard first to prevent multiple declarations
#ifndef FOCAL_H
#define FOCAL_H

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define ABSORBER_DAMPING(eco,damp) (1.-eco*damp*damp)

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
        z2waist,
        Y_at_X1, k0Ln_at_X1, theta_at_X1;
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

int set_densityInAbsorber_v2( gridConfiguration *gridCfg,
                              char absorber[], 
                              double n_e[gridCfg->Nx/2][gridCfg->Ny/2][gridCfg->Nz/2] );

int apply_absorber( gridConfiguration *gridCfg, 
                    double eco, 
                    double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] );

int apply_absorber_ref( gridConfiguration *gridCfg, 
                        double eco, 
                        double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref] );

int apply_absorber_v2( size_t N_x, size_t N_y, size_t N_z, int d_absorb, double eco, 
                       char absorber[],
                       double EB_WAVE[N_x][N_y][N_z] );

int apply_numerical_viscosity( gridConfiguration *gridCfg,
                               double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz] );

int abc_Mur_saveOldE_xdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[8][gridCfg->Ny][gridCfg->Nz] );
int abc_Mur_saveOldE_ydir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[gridCfg->Nx][8][gridCfg->Nz] );    // was [Ny,8,Nz] until 2022-01-22 (which is wrong)
int abc_Mur_saveOldE_zdir( gridConfiguration *gridCfg, 
                           double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                           double E_old[gridCfg->Nx][gridCfg->Ny][8] );

int abc_Mur_saveOldEref_xdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[8][gridCfg->Ny][gridCfg->Nz] );
int abc_Mur_saveOldEref_ydir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[gridCfg->Nx][8][gridCfg->Nz] );
int abc_Mur_saveOldEref_zdir( gridConfiguration *gridCfg, 
                              double EB_WAVE_ref[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                              double E_old[gridCfg->Nx][gridCfg->Ny][8] );

int abc_Mur_1st( gridConfiguration *gridCfg, 
                 char absorber[],
                 double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz], 
                 double E_old_xdir[8][gridCfg->Ny][gridCfg->Nz], 
                 double E_old_ydir[gridCfg->Nx][8][gridCfg->Nz], 
                 double E_old_zdir[gridCfg->Nx][gridCfg->Ny][8] );
int abc_Mur_1st_ref( gridConfiguration *gridCfg, 
                     double EB_WAVE[gridCfg->Nx][gridCfg->Ny][gridCfg->Nz_ref], 
                     double E_old_xdir[8][gridCfg->Ny][gridCfg->Nz_ref], 
                     double E_old_ydir[gridCfg->Nx][8][gridCfg->Nz_ref], 
                     double E_old_zdir[gridCfg->Nx][gridCfg->Ny][8] );

int set2zero_1D( size_t N_x, double arr_1D[N_x] );
int set2zero_3D( size_t N_x, size_t N_y, size_t N_z, double arr_3D[N_x][N_y][N_z] );

#endif  // FOCAL_H

