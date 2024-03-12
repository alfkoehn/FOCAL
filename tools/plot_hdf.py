# coding=utf-8

__author__      = 'Alf Köhn-Seemann'
__email__       = 'koehn@igvp.uni-stuttgart.de'
__copyright__   = 'University of Stuttgart'


# import standard modules
import argparse
import h5py
import matplotlib.pyplot as plt
import numpy as np
import os.path

# import modules for visualization of 3D data
from mayavi import mlab
import pyvista as py


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
    else:
        print( 'ERROR: dataset <{0}> does not exists in file <{1}>'.format( dSet_name, fname ) )
        return err_value

    h5f.close()

    return data_in
    #;}}}


def plot_simple( fname_in, dSet_name='',
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

    colScale    = 'lin'
    if colScale == 'lin':
        contLevels  = np.linspace( 0, np.amax(data2plot), 20 )[1:].tolist()
    elif colScale == 'log':
        contLevels  = np.logspace( np.log10(1e-2), np.log10(np.amax(E_abs)), 8)[3:].tolist()

    if not silent:
        print( funcName )
        print( "  status: contour levels = ", contLevels )

    fig1    = mlab.figure( bgcolor=(1,1,1),
                           fgcolor=(0,0,0), # color of axes, orientation axes, label
                           size=(800,600),
                         )

    cont_Eabs   = mlab.contour3d( data2plot, 
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


def main():
    #{{{

    # initialize parser for command line options
    parser  = argparse.ArgumentParser()
    # add optional arguments
    parser.add_argument( "-f", "--filename", type=str, default="fileout.h5",
                         help="Filename of hdf5 output file from FOCAL." )
    parser.add_argument( "-d", "--dSet_name", type=str, default="",
                         help="Dataset name to be plotted." )
    # read all argments from command line
    args        = parser.parse_args()
    fname       = args.filename
    dSet_name   = args.dSet_name

    print( "  Following arguments are set (via command line options): " )
    print( "    fname = {0}".format(fname) )
    print( "    dSet_name = {0}".format(dSet_name) )

    plot_simple(fname, dSet_name=dSet_name, silent=False)

    #}}}


if __name__ == '__main__':
    main()

