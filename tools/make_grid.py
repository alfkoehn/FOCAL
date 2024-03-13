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


def write2hdf5( params, data2save, 
                fname='grid.h5', dSet_name='dataSet',
                write_config=False
              ):
#{{{
    """
    Write into hdf5-file.

    Parameters
    ----------
    params:
    data2save:
    fname: str
        filename of hdf5-file including full path
    dSet_name: str
        name of dataset to be read from hdf5-file
    write_config:

    Returns
    -------
    """

    print( 'write2hdf5' )
    print( '    starting now to save full-wave grids into hdf5-file' )

    # create file object (append if exists, create otherwise)
    h5f = h5py.File( fname, 'a' )
    
    # optionally include configurational parameters
    if write_config:
        h5f.create_dataset( 'config/period', data=params['period'] )
        h5f.create_dataset( 'config/f_0',    data=params['f_0'] )
        h5f.create_dataset( 'config/N_z',    data=params['N_z'] )
        h5f.create_dataset( 'config/N_y',    data=params['N_y'] )

    # write array provided in function call into file
    # check if it exists
    if dSet_name in h5f.keys():
        print( '    ERROR: dataset <{0}> already exists in file, will NOT be overwritten'.format(dSet_name) )
    else:
        h5f.create_dataset( dSet_name, data=data2save, 
                            shuffle=True, compression="gzip", compression_opts=9
                          )

    h5f.close()

    print( '    done, dataset <{0}> successfully written into file {1}'.format( dSet_name, fname ) )
    #}}}
 

def make_ne_profile( ne_profile, Nx=100, Ny=70, Nz=40, 
                     fname='grid.h5', dSet_name='n_e',
                   ):
    #{{{

    # create empty array for n_e profile
    arr = np.empty( [Nx, Ny, Nz] )

    if ne_profile == 1:
        for ii in range(Nx):
            for jj in range(Ny):
                for kk in range(Nz):
                    arr[ii, jj, kk] = ii+jj+kk

    elif ne_profile == 2:
        arr[ :int(Nx/2), :, : ]  = 0
        arr[ int(Nx/2):, :, : ]  = 5

    return arr
    #}}}


def main():
    #{{{

    n_e = make_ne_profile( 2, Nx=400, Ny=300, Nz=200 )
    
    write2hdf5( [], n_e, fname='grid.h5', dSet_name='n_e' )

    #}}}


if __name__ == '__main__':
    main()

