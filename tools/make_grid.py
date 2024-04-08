# coding=utf-8

__author__      = 'Alf Köhn-Seemann'
__email__       = 'koehn@igvp.uni-stuttgart.de'
__copyright__   = 'University of Stuttgart'


# import standard modules
import argparse
import h5py
import numpy as np
import os.path

# import modules for visualization of 3D data
from mayavi import mlab


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
    data2save: numpy array
    fname: str
        filename of hdf5-file including full path
    dSet_name: str
        name of dataset to be read from hdf5-file
    write_config: bool
        include configurational full-wave parameters if true

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
        h5f.create_dataset( 'config/N_x',    data=params['N_x'] )

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


def Rx(alpha, angle='deg'):
    #{{{
    """
    Returns the rotation matrix for rotation around the x-axis. A rotation 
    of a vector would be performed by vec_new=Rx*vec_old.
    """
    if angle == 'deg':
        alpha   = np.radians(alpha)
    return np.matrix([[ 1, 0            , 0             ],
                      [ 0, np.cos(alpha), -np.sin(alpha)],
                      [ 0, np.sin(alpha), np.cos(alpha) ]])
    #}}}


def Ry(beta, angle='deg'):
    #{{{
    """
    Returns the rotation matrix for rotation around the y-axis. A rotation 
    of a vector would be performed by vec_new=Ry*vec_old.
    """
    if angle == 'deg':
        beta   = np.radians(beta)
    return np.matrix([[  np.cos(beta) , 0. , np.sin(beta ) ],
                      [  0.           , 1. , 0.            ],
                      [ -np.sin(beta) , 0. , np.cos(beta) ]])
    #}}}


def Rz(gamma, angle='deg'):
    #{{{
    """
    Returns the rotation matrix for rotation around the z-axis. A rotation 
    of a vector would be performed by vec_new=Rz*vec_old.
    """
    if angle == 'deg':
        gamma   = np.radians(gamma)
    return np.matrix([[ np.cos(gamma) , -np.sin(gamma) , 0. ],
                      [ np.sin(gamma) ,  np.cos(gamma) , 0. ],
                      [ 0             , 0.             , 1. ]])
    #}}}


def rotate_vec_via_skewing( coords_in, angle_in_degrees, 
                            rot_axis='x', rot_center=np.array([0,0,0]) ):
    #{{{
    """
    Rotate a 3-element vector via shearing/skewing.

    Parameters
    ----------
    coords_in: numpy array with 3 elements
        coordinates of vector to be rotated
    angle_in_degrees: float
        rotation angle in degrees
    rot_axis: str
        axis around which the vector
    rot_center: numpy array with 3 elements
        center around which rotation will be performed

    Returns
    -------
    numpy array with 3 elements
    """

    alpha   = np.radians(angle_in_degrees)
    
    # transform coordinate system such that rotation center is at origin
    coords_in[0] -= rot_center[0]
    coords_in[1] -= rot_center[1]
    coords_in[2] -= rot_center[2]

    # translations
    skew_1  = np.tan(alpha/2.)
    skew_2  = -np.sin(alpha)
    skew_3  = skew_1

    if rot_axis == 'x':
        x_new   = coords_in[0]
        # apply the translations (i.e. skew the vector)
        z_new   = coords_in[2] + round(skew_1*coords_in[1])
        y_new   = coords_in[1] + round(skew_2*z_new)
        z_new   = z_new + round(skew_3*y_new)

    elif rot_axis == 'y':
        y_new   = coords_in[1]
        # apply the translations (i.e. skew the vector)
        z_new   = coords_in[2] + round(skew_1*coords_in[0])
        x_new   = coords_in[0] + round(skew_2*z_new)
        z_new   = z_new + round(skew_3*x_new)

    elif rot_axis == 'z':
        z_new   = coords_in[2]
        # apply the translations (i.e. skew the vector)
        y_new   = coords_in[1] + round(skew_1*coords_in[0])
        x_new   = coords_in[0] + round(skew_2*y_new)
        y_new   = y_new + round(skew_3*x_new)

    # transform coordinate system back to original system
    x_new += rot_center[0]
    y_new += rot_center[1]
    z_new += rot_center[2]
    
    return np.array( [x_new, y_new, z_new] )
    #}}}


def rotate_arr3D_via_shearing( arr_in, angle_in_degrees, 
                               rot_axis='z', rot_center=np.array([0,0,0]) ):
    #{{{

    """
    Rotate a 3D array via shearing/skewing.

    Parameters
    ----------
    arr_in: numpy array with Nx*Ny*Nz elements
        array to be rotated
    angle_in_degrees: float
        rotation angle in degrees
    rot_axis: str
        axis around which the vector
    rot_center: numpy array with 3 elements
        center around which rotation will be performed

    Returns
    -------
    numpy array with Nx*Ny*Nz elements
    """

    # copy of original array for rotated array
    Nx, Ny, Nz  = arr_in.shape[0], arr_in.shape[1], arr_in.shape[2]
    arr_rotated = np.zeros( (Nx,Ny,Nz) )

    # note that this can probably be done more efficiently w/o for-loops
    for xx in range(Nx):
        for yy in range(Ny):
            for zz in range(Nz):
                # rotate via shearing
                coords_new = rotate_vec_via_skewing(  np.array([xx,yy,zz]), angle_in_degrees, 
                                                      rot_axis=rot_axis, 
                                                      rot_center=rot_center )
                # handle coordinates which are out of bounds after rotation
                if check_boundaries(coords_new, np.array([0,0,0]), np.array([Nx, Ny, Nz])):
                    arr_rotated[ coords_new[0], 
                                 coords_new[1], 
                                 coords_new[2] ] = arr_in[xx,yy,zz]
                
    return arr_rotated
    #}}}


def rotate_arr3D_via_rotMatrix( arr_in, angle_in_degrees,
                                rot_axis='z', rot_center=np.array([0,0,0])):
    #{{{
    """
    Rotate a 3D array via rotation matrices.

    Parameters
    ----------
    arr_in: numpy array with Nx*Ny*Nz elements
        array to be rotated
    angle_in_degrees: float
        rotation angle in degrees
    rot_axis: str
        axis around which the vector
    rot_center: numpy array with 3 elements
        center around which rotation will be performed

    Returns
    -------
    numpy array with Nx*Ny*Nz elements
    """

    # copy of original array for rotated array
    Nx, Ny, Nz  = arr_in.shape[0], arr_in.shape[1], arr_in.shape[2]
    arr_rotated = np.zeros( (Nx,Ny,Nz) )

    for xx in range(Nx):
        for yy in range(Ny):
            for zz in range(Nz):
                if arr_in[xx,yy,zz] > 0:
                    # rotate via rotation matrix
                    if rot_axis == 'x':
                        coords_new  = Rx(angle_in_degrees,angle='deg')*np.array([[xx],[yy],[zz]])
                    elif rot_axis == 'y':
                        coords_new  = Ry(angle_in_degrees,angle='deg')*np.array([[xx],[yy],[zz]])
                    elif rot_axis == 'z':
                        coords_new  = Rz(angle_in_degrees,angle='deg')*np.array([[xx],[yy],[zz]])
                    # round new coordinates to integers (as they are array indices)
                    coords_new = np.array([round(coords_new[0,0]), 
                                           round(coords_new[1,0]),
                                           round(coords_new[2,0])])
                    # check if rotated coordinates are out of boundaries
                    if check_boundaries(coords_new, 
                                        np.array([0,0,0]), np.array([Nx,Ny,Nz])):
                        arr_rotated[ coords_new[0], 
                                     coords_new[1], 
                                     coords_new[2] ] = arr_in[xx,yy,zz]
                
    return arr_rotated
    #}}}


def check_boundaries( coords_new, coords_min, coords_max ):
    #{{{

    # check if coordinates are out of bounds, to be used after rotation/translation
    # return false if out of bounds, true if ok

    # check boundaries
    if coords_new[0] < coords_min[0]    : return False#np.nan#coords_new[0] = 0
    elif coords_new[0] >= coords_max[0] : return False#np.nan#coords_new[0] = coords_max[0]-1

    if coords_new[1] < coords_min[1]    : return False#np.nan#coords_new[1] = 0
    elif coords_new[1] >= coords_max[1] : return False#np.nan#coords_new[1] = coords_max[1]-1

    if coords_new[2] < coords_min[2]    : return False#np.nan#coords_new[2] = 0
    elif coords_new[2] >= coords_max[2] : return False#np.nan#coords_new[2] = coords_max[2]-1

    return True
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
        # cuboid
        xc, yc, zc  = Nx/2, Ny/2, Nz/2
        dx, dy, dz  = Nx/4, Ny/4, Nz/4
        ne_max      = 5
        arr[:,:,:]  = 0
        arr[ round(xc-dx/2):round(xc+dx/2), 
             round(yc-dy/2):round(yc+dy/2),   
             round(zc-dz/2):round(zc+dz/2) ] = ne_max

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

    elif ne_profile == 9:
        # cuboid rotated around x-axis
        xc, yc, zc  = Nx/2, Ny/2, Nz/2
        dx, dy, dz  = Nx/4, Ny/4, Nz/4
        ne_max      = 5
        gamma       = 20    # rotation angle in degrees
        # make the original cuboid
        arr[:,:,:]  = 0
        arr[ round(xc-dx/2):round(xc+dx/2), 
             round(yc-dy/2):round(yc+dy/2),   
             round(zc-dz/2):round(zc+dz/2) ] = ne_max

        # rotate cuboid
        arr_rotated = rotate_arr3D_via_rotMatrix( arr, gamma,
                                                  rot_axis='x', 
                                                  rot_center=np.array([0,0,0]))

        # overlay original and rotated cuboid
        arr = arr + arr_rotated

    elif ne_profile == 10:
        # cuboid rotated around y-axis
        xc, yc, zc  = Nx/2, Ny/2, Nz/2
        dx, dy, dz  = Nx/4, Ny/4, Nz/4
        ne_max      = 5
        gamma       = 20    # rotation angle in degrees
        # make the original cuboid
        arr[:,:,:]  = 0
        arr[ round(xc-dx/2):round(xc+dx/2), 
             round(yc-dy/2):round(yc+dy/2),   
             round(zc-dz/2):round(zc+dz/2) ] = ne_max

        # rotate cuboid
        arr_rotated = rotate_arr3D_via_rotMatrix( arr, gamma,
                                                  rot_axis='y', 
                                                  rot_center=np.array([0,0,0]))

        # overlay original and rotated cuboid
        arr = arr + arr_rotated

    elif ne_profile == 11:
        # cuboid rotated around z-axis
        xc, yc, zc  = Nx/2, Ny/2, Nz/2
        dx, dy, dz  = Nx/4, Ny/4, Nz/4
        ne_max      = 5
        gamma       = 20    # rotation angle in degrees
        # make the original cuboid
        arr[:,:,:]  = 0
        arr[ round(xc-dx/2):round(xc+dx/2), 
             round(yc-dy/2):round(yc+dy/2),   
             round(zc-dz/2):round(zc+dz/2) ] = ne_max

        # rotate cuboid
        arr_rotated = rotate_arr3D_via_rotMatrix( arr, gamma,
                                                  rot_axis='z', 
                                                  rot_center=np.array([0,0,0]))

        # overlay original and rotated cuboid
        arr = arr + arr_rotated

    elif ne_profile == 12:
        # cuboid rotated around x-axis, with rotation realized by 3 skews
        # read here: 
        xc, yc, zc  = Nx/2, Ny/2, Nz/2
        dx, dy, dz  = Nx/4, Ny/4, Nz/4
        ne_max      = 5
        alpha       = 20    # rotation angle in degrees
        #make the original cuboid
        arr[:,:,:]  = 0
        arr[ round(xc-dx/2):round(xc+dx/2), 
             round(yc-dy/2):round(yc+dy/2),   
             round(zc-dz/2):round(zc+dz/2) ] = ne_max
        
        # rotate via shearing
        arr = rotate_arr3D_via_shearing( arr, alpha, 
                                         rot_axis='z', rot_center=np.array([0,0,0]) )


    return arr
    #}}}


def main():
    #{{{

    # initialize parser for command line options
    parser  = argparse.ArgumentParser()
    # add optional arguments
    parser.add_argument( "-f", "--filename", type=str, default="grid.h5",
                         help="Filename of hdf5 output file for FOCAL." )
    parser.add_argument( "-n", "--ne_profile", type=int, default=12,
                         help="ne_profile (look into source code for all possible options)." )

    # read all argments from command line
    args        = parser.parse_args()
    fname       = args.filename
    ne_profile  = args.ne_profile

    n_e = make_ne_profile( ne_profile, Nx=int(400/2), Ny=int(300/2), Nz=int(200/2) )
    
    write2hdf5( [], n_e, fname=fname, dSet_name='n_e' )

    #}}}


if __name__ == '__main__':
    main()

