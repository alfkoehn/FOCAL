#include "antenna.h"

static double **antField_xy = NULL;
static double **antPhaseTerms = NULL;

void init_antennaInjection( gridConfiguration *gridCfg, 
                            beamAntennaConfiguration *beamCfg ){

    //initializevalues for antenna injection
    T_wave      = 0;
    omega_t     = .0;

    //Allocate memory for AntField_xy and antPhaseTerms
    // antenna: envelope of injected field
    antField_xy = allocate2DArray(NX/2, NY/2);
    // antenna: phase terms
    antPhaseTerms = allocate2DArray(NX/2, NY/2);

    /*Initialize antenna beam*/
    printf( "starting to define antenna field...\n" );
    make_antenna_profile( gridCfg, beamCfg );
    printf( "...done defining antenna field\n" );

}

void control_antennaInjection(  gridConfiguration *gridCfg, 
                                beamAntennaConfiguration *beamCfg,
                                int t_int,
                                double EB_WAVE[NX][NY][NZ],
                                double EB_WAVE_ref[NX][NY][NZ_REF] ){

    omega_t += 2.*M_PI/period;

    // to avoid precision problems when a lot of pi's are summed up        
    if (omega_t >= 2.*M_PI) {
        omega_t    += -2.*M_PI;
        T_wave     += 1;
        //printf("status: number of oscillation periods: %d (t_int= %d) \n",T_wave,t_int);
    }

    // add source
    add_source( gridCfg, beamCfg,
                t_int,  
                EB_WAVE );
    add_source_ref( gridCfg, beamCfg,
                    t_int,  
                    EB_WAVE_ref );

}

int make_antenna_profile(   gridConfiguration *gridCfg, 
                            beamAntennaConfiguration *beamCfg ) {
//{{{
// like make_antenna_profile_3 but with previously missing optional for z2waist
// i.e. allowing now for converging beams with their waist not in the antenna plane

    int
        ii,jj;
    double
        antBeam_z_x, antBeam_z_y,
        antBeam_r_x, antBeam_r_y,
        antBeam_wx, antBeam_wy,

        antPhase_x, antPhase_y,
        antPhaseCurve_xR, antPhaseCurve_yR,
        antPhaseCurve_x, antPhaseCurve_y,
        antPhaseGouy_x, antPhaseGouy_y;

    for (ii=0 ; ii<(NX/2) ; ++ii) {
        // beam coordinate system
        antBeam_r_x  = ((double)ii-(double)ant_x/2.) * cos(antAngle_zx/180.*M_PI);
        antBeam_z_x  = ((double)ii-(double)ant_x/2.) * sin(antAngle_zx/180.*M_PI) * cos(antAngle_zy/180.*M_PI) + z2waist/2;

        // account for tilted Gauss beam
        // w(z)=w0*sqrt(1+(lambda*z/pi*w0^2)^2)
        antBeam_wx  = ant_w0x*(period/2.) * sqrt( 1. + pow( (period/2)*antBeam_z_x/( M_PI*pow(ant_w0x*(period/2), 2) ) , 2)  );

        // phase variation along beam in atenna plane
        antPhase_x  = antBeam_z_x * 2.*M_PI/(period/2.);

        // phase variation due to curvature of phase fronts
        // radius of curvature of phasefronts: R(z)=z+1/z*(pi*w0^2/lambda)^2
        antPhaseCurve_xR    = antBeam_z_x + 1./(antBeam_z_x + 1e-5) 
                                           *pow( M_PI * pow(ant_w0x*period/2., 2) / (period/2) , 2 );
        antPhaseCurve_x     = pow(antBeam_r_x,2) / (2.*antPhaseCurve_xR) * 2.*M_PI/(period/2);

        for (jj=0 ; jj<(NY/2) ; ++jj) {
            // beam coordinate system
            antBeam_r_y  = ((double)jj-(double)ant_y/2.) * cos(antAngle_zy/180.*M_PI);
            antBeam_z_y  = ((double)jj-(double)ant_y/2.) * sin(antAngle_zy/180.*M_PI) * cos(antAngle_zx/180.*M_PI) + z2waist/2;
        
            // account for tilted Gauss beam
            // w(z)=w0*sqrt(1+(lambda*z/pi*w0^2)^2)
            antBeam_wy  = ant_w0y*(period/2.) * sqrt( 1. + pow( (period/2.)*antBeam_z_y/( M_PI*pow(ant_w0y*(period/2.), 2) ) , 2)  );

            // envelope of antenna field
            antField_xy[ii][jj] = exp( -1.*pow(antBeam_r_x/antBeam_wx, 2) ) 
                                 *exp( -1.*pow(antBeam_r_y/antBeam_wy, 2) );
            // factor: w0/w(z)
            antField_xy[ii][jj] *= ant_w0x*(period/2)/antBeam_wx * ant_w0y*(period/2)/antBeam_wy;

            // phase variation along beam in atenna plane
            antPhase_y          = antBeam_z_y * 2.*M_PI/(period/2.);

            // phase variation due to curvature of phase fronts
            // radius of curvature of phasefronts: R(z)=z+1/z*(pi*w0^2/lambda)^2
            antPhaseCurve_yR    = antBeam_z_y + 1./(antBeam_z_y + 1e-5) 
                                               *pow( M_PI * pow(ant_w0y*period/2., 2) / (period/2.) , 2 );
            antPhaseCurve_y     = pow(antBeam_r_y,2) / (2.*antPhaseCurve_yR) * 2.*M_PI/(period/2.);

            // account for the Gouy-phase
            // phase_Gouy = arctan(z/z_R) 
            // with z_R = pi*w_0^2/lambda the Rayleigh range
            antPhaseGouy_x  = atan( period/2.*antBeam_z_x / (M_PI * pow(ant_w0x*period/2., 2) ) );
            antPhaseGouy_y  = atan( period/2.*antBeam_z_y / (M_PI * pow(ant_w0y*period/2., 2) ) );

                //ant_phase   = .0; <<--- extra phase-term

            antPhaseTerms[ii][jj]   = -antPhase_x 
                                      -antPhase_y
                                      -antPhaseCurve_x
                                      -antPhaseCurve_y
                                      -antPhaseGouy_x
                                      -antPhaseGouy_y;
        }
    }

    return EXIT_SUCCESS;
}//}}}


int add_source( gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg, 
                int t_int,
                double EB_WAVE[NX][NY][NZ] ) {
//{{{

    size_t
        ii, jj;
    double
        fact1_Hansen_corr,
        fact2_Hansen_corr,
        theta_rad,
        t_rise, 
        source;

    // slowly increase field in time 
    t_rise  = antenna_field_rampup( rampUpMethod, period, t_int );

    if ( exc_signal == 1 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii, jj, source)
        for ( ii=2 ; ii<NX ; ii+=2 ) {
            for ( jj=2 ; jj<NY ; jj+=2 ) {
                // note: for X-mode injection, switch cos and sin of source_1 and source_2
                //source      = sin(omega_t - aux - curve + GouyPhase_beam + ant_phase/180.*M_PI ) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                // Ex
                EB_WAVE[ii+1][jj  ][ant_z]   += source;
            }
        }
    } else if ( exc_signal == 2) {
#pragma omp parallel for collapse(2) default(shared) private(ii, jj, source)
        for ( ii=2 ; ii<NX ; ii+=2 ) {
            for ( jj=2 ; jj<NY ; jj+=2 ) {
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                // Bx
                EB_WAVE[ii  ][jj+1][ant_z+1]   += source;
            }
        }
    } else if ( exc_signal == 3) {
#pragma omp parallel for collapse(2) default(shared) private(ii, jj, source)
        for ( ii=2 ; ii<NX ; ii+=2 ) {
            for ( jj=2 ; jj<NY ; jj+=2 ) {
                // note: for X-mode injection, switch cos and sin of source
                //       or, add/subtract pi/2 in sine for Bx 
                // Ex
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                EB_WAVE[ii+1][jj  ][ant_z]   += source;
                // Bx
                //source  = cos(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)] + M_PI/2.) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                EB_WAVE[ii  ][jj+1][ant_z+1] += source*(1.41)*0.;
            }
        }
    } else if ( exc_signal == 4) {
        // elliptically polarized for optimum O-SX conversion using Hansen's
        // formula for calculating ratio of wave electric fields 

        // calculate factor defining ratio of perpendicular E-fields according to Hansen
        //theta_rad           = beamCfg->antAngle_zx/180. * M_PI;
        theta_rad           = theta_at_X1/180. * M_PI;
        fact1_Hansen_corr   = antenna_calcHansenExEy_O( theta_rad, Y_at_X1 );
        fact2_Hansen_corr   = -1.*cos(theta_rad)/sin(theta_rad) * .0;

        if (t_int < 1) {
            printf( "|E_x/E_y| = %f\n", fact1_Hansen_corr );
            printf( " E_x/E_z  = %f\n", fact2_Hansen_corr );
        }

#pragma omp parallel for collapse(2) default(shared) private(ii, jj, source)
        for ( ii=2 ; ii<NX ; ii+=2 ) {
            for ( jj=2 ; jj<NY ; jj+=2 ) {
                // Ex
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                EB_WAVE[ii+1][jj  ][ant_z  ] += source;
                // Ey
                //source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)] + M_PI/2.) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                source  = cos(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                EB_WAVE[ii  ][jj+1][ant_z  ] += source/fact1_Hansen_corr;
                // Ez
                //source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                //EB_WAVE[ii  ][jj  ][beamCfg->ant_z+1] += source/fact2_Hansen_corr;
                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                // WARNING: Ez has been switched off 
                //          ==> no O-mode injection for angled injection possible this way
                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            }
        }
    } else if ( exc_signal == 5) {
        // linearly polarized beam for oblique injection (added: 2023-02-13)
#pragma omp parallel for collapse(2) default(shared) private(ii, jj, source)
        for ( ii=2 ; ii<NX ; ii+=2 ) {
            for ( jj=2 ; jj<NY ; jj+=2 ) {
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                // Ex
                EB_WAVE[ii+1][jj  ][ant_z  ] += source * (1.*cos(antAngle_zx/180.*M_PI));
            }
        }
    }

    return EXIT_SUCCESS;
}//}}}


int add_source_ref( gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg, 
                    int t_int,  
                    double EB_WAVE[NX][NY][NZ_REF] ) {
//{{{

    size_t
        ii, jj;
    double
        fact1_Hansen_corr,
        fact2_Hansen_corr,
        theta_rad,
        t_rise, 
        source;

    // slowly increase field in time 
    t_rise  = antenna_field_rampup( rampUpMethod, period, t_int );

    if ( exc_signal == 1 ) {
#pragma omp parallel for collapse(2) default(shared) private(ii, jj, source)
        for ( ii=2 ; ii<NX ; ii+=2 ) {
            for ( jj=2 ; jj<NY ; jj+=2 ) {
                // note: for X-mode injection, switch cos and sin of source_1 and source_2
                //source      = sin(omega_t - aux - curve + GouyPhase_beam + ant_phase/180.*M_PI ) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                // Ex
                EB_WAVE[ii+1][jj  ][ant_z  ]   += source;
            }
        }
    } else if ( exc_signal == 2) {
        t_rise  = 1. - exp( -1*pow( ((double)(t_int)/period), 2 )/100. );
#pragma omp parallel for collapse(2) default(shared) private(ii, jj, source)
        for ( ii=2 ; ii<NX ; ii+=2 ) {
            for ( jj=2 ; jj<NY ; jj+=2 ) {
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                // Bx
                EB_WAVE[ii  ][jj+1][ant_z+1]   += source;
            }
        }
    } else if ( exc_signal == 3) {
#pragma omp parallel for collapse(2) default(shared) private(ii, jj, source)
        for ( ii=2 ; ii<NX ; ii+=2 ) {
            for ( jj=2 ; jj<NY ; jj+=2 ) {
                // note: for X-mode injection, switch cos and sin of source_1 and source_2
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                // Ex
                EB_WAVE[ii+1][jj  ][ant_z  ]   += source;
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)] + M_PI/2.) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                // Bx
                EB_WAVE[ii  ][jj+1][ant_z+1]   += source*(1.41);
            }
        }
    } else if ( exc_signal == 4) {
        // elliptically polarized for optimum O-SX conversion using Hansen's
        // formula for calculating ratio of wave electric fields 

        //theta_rad           = beamCfg->antAngle_zx/180. * M_PI;
        theta_rad           = theta_at_X1/180. * M_PI;
        fact1_Hansen_corr   = antenna_calcHansenExEy_O( theta_rad, Y_at_X1 );
        fact2_Hansen_corr   = -1.*cos(theta_rad)/sin(theta_rad) * .0;

        if (t_int < 1) {
            printf( "|E_x/E_y| = %f\n", fact1_Hansen_corr );
            printf( " E_x/E_z  = %f\n", fact2_Hansen_corr );
        }

#pragma omp parallel for collapse(2) default(shared) private(ii, jj, source)
        for ( ii=2 ; ii<NX ; ii+=2 ) {
            for ( jj=2 ; jj<NY ; jj+=2 ) {
                // Ex
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                EB_WAVE[ii+1][jj  ][ant_z  ] += source;
                // Ey
                //source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)] + M_PI/2.) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                source  = cos(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                EB_WAVE[ii  ][jj+1][ant_z  ] += source/fact1_Hansen_corr;
                // Ez
                //source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                //EB_WAVE[ii  ][jj  ][beamCfg->ant_z+1] += source/fact2_Hansen_corr;
                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                // WARNING: Ez has been switched off 
                //          ==> no O-mode injection for angled injection possible this way
                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            }
        }
    } else if ( exc_signal == 5) {
        // linearly polarized beam for oblique injection (added: 2023-02-13)
#pragma omp parallel for collapse(2) default(shared) private(ii, jj, source)
        for ( ii=2 ; ii<NX ; ii+=2 ) {
            for ( jj=2 ; jj<NY ; jj+=2 ) {
                source  = sin(omega_t + antPhaseTerms[(ii/2)][(jj/2)]) * t_rise * antField_xy[(ii/2)][(jj/2)] ;
                // Ex
                EB_WAVE[ii+1][jj  ][ant_z  ] += source * (1.*cos(antAngle_zx/180.*M_PI));
            }
        }
    }


    return EXIT_SUCCESS;
}//}}}


double antenna_field_rampup( int RampUpMethod, double Period, int t_int ){
//{{{

    // If the amplitude of the wave electric field is increased too fast,
    // higher harmonics can in principle be excited which can result in an 
    // oscillating behaviour in the detected power when analyzing mode
    // conversion scenarios. To avoid this, the field is increased slowly
    // in time. To achieve this, this function returns a factor which is to
    // be multiplied with the antenna field added to the full-wave grid (i.e.
    // the antenna) at each time step. This factor, called t_rise, will be 0
    // at the beginning and ramped-up to 1 exponentially. The speed of the 
    // ramp-up can be varied with the parameter tau, for tau=100, this 
    // corresponds to needed roughly 30 oscillation periods to reach one. 
    
    double
        t_rise,
        tau;

    tau = 100.;

    if (RampUpMethod == 1) {
        // exponential increase reaching 1 after roughly 30 oscillation periods
        t_rise  = 1. - exp( -1*pow( ((double)(t_int)/Period), 2 )/tau );
    } else {
        printf( "antenna_field_rampup: WARNING, rampUpMethod %d does not exist\n", RampUpMethod );
        printf( "                      ==> no smooth ramp-up, just set instantly to 1\n" );
        t_rise  = 1.;
    }

    return t_rise;
} //}}}


double antenna_calcHansenExEy_O( double theta_rad, double Y ){
    //{{{
    //
    // TODO: requires tests if it actually works
    //       compare this value with value deduced from python scripts

    double
        ExEy;

    ExEy    = .5*(Y*pow(sin(theta_rad),2) 
    //ExEy    = .5*(Y*Y*pow(sin(theta_rad),2) // <--- is this correct ? TODO: check! also check sign!!
                  +sqrt( Y*Y*pow(sin(theta_rad),4) + 4*pow(cos(theta_rad),2) )
                 );

    ExEy    = 0.52;

    return ExEy;
}//}}}
