import matplotlib.pyplot as plt
from mayavi import mlab
import numpy as np
from tvtk.api import tvtk

import matplotlib
import mayavi
import sys

# check versions
print( 'matplotlib version: ', matplotlib.__version__ )
print( 'mayavi version    : ', mayavi.__version__ )
print( 'numpy version     : ', np.__version__ )
print( 'python version    : ', sys.version_info )

# create the grid
x       = np.linspace(0,np.pi,100)
y       = np.linspace(0,2*np.pi,200)
z       = np.linspace(0,2*np.pi,400)
X, Y, Z = np.meshgrid(x, y, z)

# create the data to be plotted
dat = np.cos(X) * np.sin(Y) * np.cos(Z)

fig1    = mlab.figure( bgcolor=(1,1,1), 
                       fgcolor=(0,0,0),
                       size=(800,600), 
                       )

# make a 3D contour plot
cont3d  = mlab.contour3d( dat, 
                          transparent=True, opacity=.4,
                          figure=fig1,
                          colormap='bwr'
                        )
mlab.outline(cont3d)
mlab.show()
