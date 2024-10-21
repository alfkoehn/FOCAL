# coding=utf-8

__author__      = 'Alf Köhn-Seemann'
__email__       = 'koehn@igvp.uni-stuttgart.de'
__copyright__   = 'University of Stuttgart'


# import standard modules
import argparse
import h5py
import numpy as np
import os.path
import scipy.constants as consts

# import modules for visualization of 3D data
from mayavi import mlab


def readhdf5( fname, dSet_name ):
    #;{{{
    """
    Open hdf5-file, extract one dataset, and return it.
    Copied from little_helper.py on 2024-03-12.

    Parameters
    ----------
    fname: str
        filename of hdf5-file including full path
    dSet_name: str
        name of dataset to be read from hdf5-file

    Returns
    -------
    numpy array
    """

    err_value = 1

    # if file exists, open it for reading
    if os.path.isfile( fname ):
        h5f = h5py.File( fname, 'r' )
    else:
        print( 'ERROR: cannot read following file: {0}'.format( fname ))
        return err_value

    # trying to read dataset
    if dSet_name in h5f:
        data_in = h5f[ dSet_name ][()]
        # note that the dataset.value attribute was deprecated, i.e. 
        # data_in = h5f.get(dSet_name).value is no longer working
    else:
        print( 'ERROR: dataset <{0}> does not exists in file <{1}>'.format( dSet_name, fname ) )
        return err_value

    h5f.close()

    return data_in
    #;}}}


def calc_wpe( density ):
#;{{{
    """
    Calculate the electron plasma angular frequency.

    w_pe = sqrt( n_e * e^2 / (epsilon_0 * m_e) )

    Parameters
    ----------
    density: float
        electron plasma density in m^-3

    Returns
    -------
    float
        angular frequency in rad s^-1
    """

    return np.sqrt( density * consts.e**2 / (consts.epsilon_0 * consts.m_e) )

#;}}}


def calc_wce( B0 ):
#;{{{
    """
    Calculate the electron cyclotron angular frequency.
    
    w_ce = e*B0/m_e

    Parameters
    ----------
    B0: float
        magnetic field strength in Tesla

    Returns
    -------
    float
        angular frequency in rad s^-1
    """

    return consts.e*B0/consts.m_e
#;}}}


def calc_wci( B0=1., m_i=consts.m_p, Z=1 ):
#{{{
    """
    Calculate the ion cyclotron angular frequency.

    w_ci = e * B0/m_i

    Parameters
    ----------
    B0: float
        magnetic field strength in Tesla
    m_i: float
        ion mass in kg
    Z: int
        number of positive or negative charges of ion
        NOTE: currently not implemented

    Returns
    -------
    float
        angular frequency in rad s^-1
    """

    return consts.e*B0/m_i
#}}}


def calc_wL( B0=1., density=1e20):
#;{{{
    """
    Calculate the left-hand angular frequency w_L.

    w_L = sqrt( 0.25* w_ce^2 + w_pe^2 ) - 0.5*w_ce

    Parameters
    ----------
    B0: float
        magnetic field strength in Tesla
    density: float
        electron plasma density in m^-3

    Returns
    -------
    float
        angular frequency in rad s^-1
    """

    # calculate the characteristic frequency needed
    w_ce = calc_wce( B0 )
    w_pe = calc_wpe( density)

    return np.sqrt( .25*w_ce**2 + w_pe**2 ) - .5*w_ce
#;}}}


def calc_wR( B0=1., density=1e20):
#;{{{
    """
    Calculates the right-hand angular frequency w_R

    w_R = sqrt( 0.25* w_ce^2 + w_pe^2 ) + 0.5*w_ce

    Parameters
    ----------
    B0: float
        magnetic field strength in Tesla
    density: float
        electron plasma density in m^-3

    Returns
    -------
    float
        angular frequency in rad s^-1
    """

    # calculate the characteristic frequency needed
    w_ce = calc_wce( B0 )
    w_pe = calc_wpe( density)

    return np.sqrt( .25*w_ce**2 + w_pe**2 ) + .5*w_ce
#;}}}


def calc_wUH( B0=1., density=1e20):
#;{{{
    """
    Calculate the upper hybrid angular frequency w_UH.

    w_UH = sqrt( w_pe^2 + w_ce^2 )

    Parameters
    ----------
    B0: float
        magnetic field strength in Tesla
    density: float
        electron plasma density in m^-3

    Returns
    -------
    float
        angular frequency in rad s^-1
    """

    # calculate the characteristic frequency needed
    w_ce = calc_wce( B0 )
    w_pe = calc_wpe( density)

    return np.sqrt( w_pe**2 + w_ce**2 ) 
#;}}}


def plot_simple( fname_in, dSet_name='',
                 contLevels=20, 
                 plotReductionLevel=4,
                 fname_out='',
                 silent=True, 
               ):
    #{{{

    funcName    = 'plot_simple'

    # check keywords
    if len(dSet_name) == 0:
        print( funcName )
        print( "  ERROR: you did not set parameter 'dSet_name' " )
        print( "         will exit now" )
        return

    data2plot   = readhdf5( fname_in, dSet_name )

    print("dataset-name = {0}, min = {1}, max = {2}".format(dSet_name, np.amin(data2plot), np.amax(data2plot)) )

    colScale    = 'lin'
    if colScale == 'lin':
        contLevels  = np.linspace( 0, np.amax(data2plot), contLevels )[1:].tolist()
        #contLevels  = np.linspace( np.amin(data2plot), np.amax(data2plot), contLevels )[1:].tolist()
    elif colScale == 'log':
        contLevels  = np.logspace( np.log10(1e-2), np.log10(np.amax(E_abs)), 8)[3:].tolist()

    if not silent:
        print( funcName )
        print( "  status: contour levels = ", contLevels )

    fig1    = mlab.figure( bgcolor=(1,1,1),
                           fgcolor=(0,0,0), # color of axes, orientation axes, label
                           size=(800,600),
                         )

    cont_Eabs   = mlab.contour3d( data2plot[::plotReductionLevel,::plotReductionLevel,::plotReductionLevel], 
                                  contours=contLevels,
                                  transparent=True, opacity=.4,
                                  figure=fig1
                                )

    # create an axes instance to modify some of its properties afterwards
    ax1 = mlab.axes( nb_labels=4,
                     extent=[1, data2plot.shape[0], 
                             1, data2plot.shape[1],
                             1, data2plot.shape[2] ],
                   )
    mlab.outline(ax1)
    ax1.axes.label_format   = '%.0f'
    # labels can also be set via mlab.xlabel('x')
    ax1.axes.x_label    = 'x'
    ax1.axes.y_label    = 'y'
    ax1.axes.z_label    = 'z'

    # set initial viewing angle
    # azimuth:   angle subtended by position vector on a sphere projected onto x-y plane with the x-axis, 0...360
    # elevation: zenith angle, i.e. angle subtended by position vector and the z-axis, 0...180
    mlab.view( azimuth=290, elevation=80 )

    if len(fname_out) > 0:
        mlab.savefig( fname_out )
    else:
        mlab.show()

    #}}}


def plot_fullwave( fname_in, fname_plot='',
                   f_0=28e9,
                   t_int=0,
                   N_contLevels=20, 
                   include_absorbers=False, oplot_dens_projection=False,
                   scale_axes_to_meters=False,
                   oplot_Efieldcut=None,
                 ):
#{{{

    # set some general properties of plot
    if include_absorbers:
        plot_abs    = 'plane'
    else:
        plot_abs    = 'not'

    # read configurational data
    period      = readhdf5( fname_in, 'config/period' )[0]
    d_absorb    = readhdf5( fname_in, 'config/d_absorb' )[0]

    l_0         = consts.c/f_0

    print( 'configurational parameters: ' )
    print( '    period  : {0}'.format(period) )
    print( '    d_absorb: {0}'.format(d_absorb) )

    # scale to Yee-grid
    d_absorb    = int(round(d_absorb/2))
    print( '  scaled to physical (Yee-cell) grid:' )
    print( '    d_absorb: {0}'.format(d_absorb) )

    #Ex  = readhdf5( fname_in, 'Ex')
    #Ey  = readhdf5( fname_in, 'Ey')
    #Ez  = readhdf5( fname_in, 'Ez')
    density = readhdf5( fname_in, 'n_e')

    B0_abs  = np.sqrt( readhdf5( fname_in, 'B0x')**2 
                      +readhdf5( fname_in, 'B0y')**2 
                      +readhdf5( fname_in, 'B0z')**2 )

    #E_abs   = np.sqrt( Ex**2 + Ey**2 + Ez**2 )

    if t_int == 0:
        dSet_name   = 'E_abs__tint00161'
    else:
        dSet_name   = 'E_abs__tint{0:05d}'.format(t_int)
    E_abs   = readhdf5( fname_in, dSet_name )
    print( dSet_name )
    print( E_abs.shape )

    if isinstance(density, np.ndarray):
        Nx      = density.shape[0]
        Ny      = density.shape[1]
        Nz      = density.shape[2]
    else:
        Nx      = E_abs.shape[0]
        Ny      = E_abs.shape[1]
        Nz      = E_abs.shape[2]
    ant_x   = Nx/2 
    ant_y   = Ny/2
    ant_z   = d_absorb + 2

    print( 'min|max(n_e)   = {0}|{1}'.format(np.amin(density), np.amax(density)) )
    print( 'min|max(|B_0|) = {0}|{1}'.format(np.amin(B0_abs), np.amax(B0_abs)) )
    print( 'min|max(E_abs) = {0}|{1}'.format(np.amin(E_abs), np.amax(E_abs)) )

    # optionally, create spatial coordinate axes scaled to meters
    if scale_axes_to_meters:
        if include_absorbers:
            x0  = 0. - d_absorb/(period/2)*l_0
            y0  = 0. - d_absorb/(period/2)*l_0
            z0  = 0. - d_absorb/(period/2)*l_0
        else:
            x0  = .0
            y0  = .0
            z0  = .0
        x1  = x0 + Nx/(period/2)*l_0
        y1  = y0 + Ny/(period/2)*l_0
        z1  = z0 + Nz/(period/2)*l_0
        print( 'x0 = {0}, x1 = {1}, Nx = {2}'.format(x0, x1, Nx) )
        print( 'y0 = {0}, y1 = {1}, Ny = {2}'.format(y0, y1, Ny) )
        print( 'z0 = {0}, z1 = {1}, Nz = {2}'.format(z0, z1, Nz) )
        #xVals   = np.linspace( x0, x1, Nx )
        #yVals   = np.linspace( y0, y1, Ny )
        #zVals   = np.linspace( z0, z1, Nz )
        X, Y, Z = np.mgrid[ x0:x1:Nx*1j, y0:y1:Ny*1j, z0:z1:Nz*1j ]
    else:
        X, Y, Z = np.mgrid[ 1:Nx:Nx*1j, 1:Ny:Ny*1j, 1:Nz:Nz*1j ]

    if not include_absorbers:
        print( 'absorbers will not be included in plot, will actually be cut from data' )
        print( 'E_abs.shape = {0}, d_absorb = {1}, Nx = {2}, Ny = {3}, Nz = {4}, ant_z = {5}'.format(
            E_abs.shape, d_absorb, Nx, Ny, Nz, ant_z) )
        if isinstance(density, np.ndarray):
            density = density[d_absorb:(Nx-d_absorb), 
                              d_absorb:(Ny-d_absorb), 
                              ant_z:(Nz-d_absorb) ]
        E_abs   = E_abs[  d_absorb:(Nx-d_absorb), 
                          d_absorb:(Ny-d_absorb), 
                          ant_z:(Nz-d_absorb) ]
        X   = X[      d_absorb:(Nx-d_absorb), 
                      d_absorb:(Ny-d_absorb), 
                      ant_z:(Nz-d_absorb) ]
        Y   = Y[      d_absorb:(Nx-d_absorb), 
                      d_absorb:(Ny-d_absorb), 
                      ant_z:(Nz-d_absorb) ]
        Z   = Z[      d_absorb:(Nx-d_absorb), 
                      d_absorb:(Ny-d_absorb), 
                      ant_z:(Nz-d_absorb) ]
        if isinstance(B0_abs, np.ndarray):
            B0_abs  = B0_abs[ d_absorb:(Nx-d_absorb), 
                              d_absorb:(Ny-d_absorb), 
                              ant_z:(Nz-d_absorb) ]


    # |E| = sqrt(Ex^2 + Ey^2 + Ez^2)
    contLevels  = np.linspace( 0, np.amax(E_abs), N_contLevels )[1:].tolist()
    #contLevels  = np.linspace( 0, .1*np.amax(E_abs), 10 )[1:].tolist()
    #contLevels  = np.logspace( np.log10(1e-2), np.log10(np.amax(E_abs)), 8 )[3:].tolist()

    print( 'contour levels: ', contLevels )

    fig1    = mlab.figure( bgcolor=(1,1,1), 
                           fgcolor=(0,0,0),     # color of axes, orientation axes, labels, ...
                           size=(800, 600), 
                         )

    cont_Eabs   = mlab.contour3d( #X, Y, Z,
                                  E_abs, contours=contLevels,
                                  transparent=True, opacity=.4,
                                  figure=fig1
                                )

    if oplot_Efieldcut:
        if include_absorbers:
            slice_x1    = ant_x
            slice_z1    = ant_z
        else:
            slice_x1    = ant_x - d_absorb
            slice_z1    = ant_z - d_absorb
        if 'x1' in oplot_Efieldcut:
            slice_Eabs  = mlab.volume_slice( E_abs,
                                             slice_index=slice_x1,
                                             plane_orientation='x_axes',
                                             figure=fig1
                                           )
        if 'z1' in oplot_Efieldcut:
            slice_Eabs  = mlab.volume_slice( E_abs,
                                             slice_index=slice_z1,
                                             plane_orientation='z_axes',
                                             figure=fig1
                                           )


    if isinstance(density, np.ndarray) and (np.amin(density) != np.amax(density)):
        cont_levels = np.arange( 1, np.amax(density), 1 ).tolist()
        cont_dens   = mlab.contour3d( #X, Y, Z, 
                                      density, 
                                      contours=cont_levels, 
                                      opacity=.3, 
                                      color=(1,0,0),
                                      #colormap='gist_yarg',
                                      figure=fig1
                                    )
        # optionally, overplot projection of density profile onto y=const plane
        if oplot_dens_projection:
            # plotting onto y=const plane, x-coordinate=density, z-coordinate=z
            #dens_line   = mlab.plot3d( density[0,0,:]*30+1, 
            #                           density[0,0,:]*0+(Ny-2*d_absorb),
            # plotting onto x=const plane, y-coordinate=density, z-coordinate=z
            dens_line   = mlab.plot3d( density[0,0,:]*0+1,  
                                       #density[0,0,:]*30+1,
                                       density[int(Nx/2),int(Ny/2),:]*30+1,
                                       np.linspace(1, Nz-2*d_absorb-3, Nz-2*d_absorb-2),
#                                   color=(0,0,0),
                                       line_width=10,
                                       tube_radius=1,
                                       figure=fig1
                                     )

            # oplot a slice of the density onto the x-axis
            slice_ne    = mlab.volume_slice( density,
                                             plane_orientation='x_axes',
                                             opacity=.4,
                                             figure=fig1
                                           )

    if isinstance(B0_abs, np.ndarray) and (np.amax(B0_abs) > .0):
    #if isinstance(density, np.ndarray) and (np.amin(B0_abs) != np.amax(B0_abs)):
        B0_res  = f_0*2*np.pi / consts.e * consts.m_e
        n_Ocut  = (f_0*2*np.pi)**2 / consts.e**2 * consts.m_e * consts.epsilon_0
        # right-hand cut-off
        RH_norm = calc_wR( B0=B0_abs*B0_res, density=(density*n_Ocut) )/(f_0*2*np.pi)
        print( 'min|max(RH_norm)   = {0}|{1}'.format(np.amin(RH_norm), np.amax(RH_norm)) )
        cont_dens   = mlab.contour3d( RH_norm, 
                                      contours=[1], 
                                      opacity=.3, 
                                      color=(1,0,0),
                                      figure=fig1
                                    )
        # upper-hybrid resonance
        UH_norm = np.sqrt( density + B0_abs**2 )
        print( 'min|max(UH_norm)   = {0}|{1}'.format(np.amin(UH_norm), np.amax(UH_norm)) )
        cont_dens   = mlab.contour3d( UH_norm, 
                                      contours=[1], 
                                      opacity=.3, 
                                      color=(0,.5,0),
                                      figure=fig1
                                    )


    if plot_abs == 'plane':
        absorber_plane  = density*.0

        # x1, x2 absorber plane
        absorber_plane[d_absorb,:,:]        = 1
        absorber_plane[(Nx-d_absorb),:,:]   = 1
        absorber_plane[:,d_absorb,:]        = 1
        absorber_plane[:,(Ny-d_absorb),:]   = 1
        absorber_plane[:,:,d_absorb]        = 1
        absorber_plane[:,:,(Nz-d_absorb)]   = 1
        cont_abs_x1   = mlab.contour3d( absorber_plane, contours=[1], 
                                        opacity=.1, colormap='gist_yarg',
                                        figure=fig1
                                      )


    # create an axes instance in order to modify some of its properties afterwards
    ax1 = mlab.axes( nb_labels=4, 
                     extent=[1, E_abs.shape[0], 1, E_abs.shape[1], 1, E_abs.shape[2] ],
            )
    #mlab.outline(cont_Eabs)
    #mlab.outline(cont_dens)
    mlab.outline(ax1 )
    ax1.axes.label_format   = '%.0f'
    # labels can also be set using  mlab.xlabel('x')
    ax1.axes.x_label        = 'x'
    ax1.axes.y_label        = 'y'
    ax1.axes.z_label        = 'z'

    # set initial viewing angle
    # azimuth:   angle subtended by position vector on a sphere projected onto x-y plane with the x-axis, 0...360
    # elevation: zenith angle, i.e. angle subtended by position vector and the z-axis, 0...180
    mlab.view( azimuth=290, elevation=80 )

    if len(fname_plot) > 0:
        mlab.savefig( fname_plot )
    else:
        mlab.show()

#}}}


def main():
    #{{{

    # initialize parser for command line options
    parser  = argparse.ArgumentParser()
    # add optional arguments
    parser.add_argument( "-f", "--filename", type=str, default="fileout.h5",
                         help="Filename of hdf5 output file from FOCAL." )
    parser.add_argument( "-d", "--dSet_name", type=str, default="",
                         help="Dataset name to be plotted." )
    parser.add_argument( "-c", "--contLevels", type=int, default=20,
                         help="Number of contour levels used." )
    parser.add_argument( "-r", "--plotReductionLevel", type=int, default=4,
                         help="Reduce resolution for 3D plot by this amount." )
    parser.add_argument( "-p", "--plot_type", type=int, default=1,
                         help="Plot type." )
    parser.add_argument( "-t", "--time", type=int, default=0,
                         help="Timestep." )

    # read all argments from command line
    args                = parser.parse_args()
    fname               = args.filename
    dSet_name           = args.dSet_name
    contLevels          = args.contLevels
    plotReductionLevel  = args.plotReductionLevel
    plot_type           = args.plot_type
    t_int               = args.time

    print( "  Following arguments are set (via command line options): " )
    print( "    fname = {0}".format(fname) )
    print( "    dSet_name = {0}".format(dSet_name) )

    if plot_type == 1:
        plot_simple(fname, dSet_name=dSet_name, 
                    contLevels=contLevels, plotReductionLevel=plotReductionLevel, 
                    silent=False)
    elif plot_type == 2:
        plot_fullwave( fname, t_int=t_int, include_absorbers=False, 
                       oplot_dens_projection=False,
                       N_contLevels=contLevels,
                       #oplot_Efieldcut='x1z1',
                     )

    #}}}


if __name__ == '__main__':
    main()

