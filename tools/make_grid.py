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
        # mirror/jump parallel to yz-plane
        arr[ :int(Nx/2), :, : ]  = 0
        arr[ int(Nx/2):, :, : ]  = 5
    elif ne_profile == 3:
        # mirror/jump parallel to xz-plane
        arr[ :, :int(Ny/2), : ]  = 0
        arr[ :, int(Ny/2):, : ]  = 5
    elif ne_profile == 4:
        # mirror/jump parallel to xy-plane
        arr[ :, :, :int(Nz/2) ]  = 0
        arr[ :, :, int(Nz/2): ]  = 5

    elif ne_profile == 5:
        # tilted xy-plane (rotated around y-axis)
        # z = m*x + b = dz/dx + z0 = (z1-z0)/dx + z0
        z0  = round(Nz/2 - Nz/5)
        z1  = round(Nz/2 + Nz/5)
        dx  = Nx
        z   = lambda x : (z1-z0)/dx*x + z0
        ne_max  = 5
        for ii in range(Nx):
            for kk in range(Nz):
                if kk < z(ii):
                    arr[ii, :, round(kk)]  = 0
                else:
                    arr[ii, :, round(kk)]  = ne_max
    elif ne_profile == 6:
        # tilted xy-plane (rotated around x-axis)
        # z = m*y + b = dz/dy + z0 = (z1-z0)/dy + z0
        z0  = round(Nz/2 - Nz/5)
        z1  = round(Nz/2 + Nz/5)
        dy  = Ny
        z   = lambda y : (z1-z0)/dy*y + z0
        ne_max  = 5
        for jj in range(Ny):
            for kk in range(Nz):
                if kk < z(jj):
                    arr[:, jj, round(kk)]  = 0
                else:
                    arr[:, jj, round(kk)]  = ne_max

    elif ne_profile == 7:
        # cube
        xc  = Nx/2
        yc  = Ny/2
        zc  = Nz/2
        dx  = Nx/4
        dy  = Ny/4
        dz  = Nz/4
        ne_max      = 5
        arr[:,:,:]  = 0
        # note that this can be done more efficient (i.e. w/o a for-loop)
        for ii in range(Nx):
            for jj in range(Ny):
                for kk in range(Nz):
                    if (    (ii > (xc-dx/2) and ii < (xc+dx/2))
                        and (jj > (yc-dy/2) and jj < (yc+dy/2))
                        and (kk > (zc-dz/2) and kk < (zc+dz/2))
                       ):
                        arr[ii,jj,kk]   = ne_max

    elif ne_profile == 8:
        # sphere
        xc  = Nx/2
        yc  = Ny/2
        zc  = Nz/2
        dr  = Nz/4
        ne_max      = 5
        arr[:,:,:]  = 0
        for ii in range(Nx):
            for jj in range(Ny):
                for kk in range(Nz):
                    r = np.sqrt( (xc-ii)**2 + (yc-jj)**2 + (zc-kk)**2 )
                    if r < dr:
                        arr[ii,jj,kk]   = ne_max

    return arr
    #}}}


def main():
    #{{{

    n_e = make_ne_profile( 8, Nx=int(400/2), Ny=int(300/2), Nz=int(200/2) )
    
    write2hdf5( [], n_e, fname='grid.h5', dSet_name='n_e' )

    #}}}


if __name__ == '__main__':
    main()

