#!/usr/bin/env python
from math import pi
import numpy as np
from numpy import array
from pyhdf.SD import *
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
from matplotlib import cm
from sympy.vector import CoordSys3D,Del
from sympy import *
import os
import numpy.ma as ma

from sympy import pretty_print as pp, latex
from sympy.abc import a, b, n

from matplotlib.pyplot import figure
import matplotlib as mpl

mpl.use('TkAgg')

mpl.rcParams['xtick.major.pad'] = 6
mpl.rcParams['ytick.major.pad'] = 6

num_of_pics = -1

fig, axs = plt.subplots(1, 3, sharey=True)
fig.subplots_adjust(wspace=0)

axs = axs.ravel()

for filename in os.listdir("./"):
    if filename.startswith("hdfra."):

        num_of_pics = num_of_pics + 1

        fo = SD(filename)
        fo_attr=fo.attributes()
        fo_dsets=fo.datasets()

        ds0=fo.select(0)

        ds1=fo.select(1)


        ds2=fo.select(2)


        ds7=fo.select('Data-Set-2')
        ds8=fo.select('Data-Set-3')
        ds9=fo.select('Data-Set-4')
        ds5=fo.select('Data-Set-5')
        ds6=fo.select('Data-Set-6')        

        phi=ds0.get()
        theta=ds1.get()
        r=ds2.get()   
        den=ds5[0]
        ene=ds6[0]
        v1=ds7[0]
        v2=ds8[0]
        v3=ds9[0]

        igrid=len(theta); jgrid=len(r)
     
        dtheta = np.zeros_like(theta)
        dtheta[1:] = theta[1:] - theta[0:-1]

        M=np.zeros((igrid,jgrid),dtype=float)

        dm=np.zeros((igrid,jgrid),dtype=float)
        m=np.zeros((igrid,jgrid),dtype=float)
        for i in range(0,igrid):
            for j in range(0,jgrid):
        #rmat=np.matrix(r)
        #denmat=np.matrix(den)
        #v1mat=np.matrix(v1)
        #thetamat=np.matrix(theta)
        #dthetamat=np.matrix(dtheta)
                dm[i,j] = (2 * pi * r[j]**2) * den[i,j] * v1[i,j] * np.sin(theta[i]) * dtheta[i]
                m[i,j]= np.sum(dm[i,j],axis=0)
                m[i,j]=abs(m[i,j])
       

        m = m.reshape((120,30))

        xx = np.zeros((igrid, jgrid), dtype=float)
        yy = np.zeros((igrid, jgrid), dtype=float)
        for i in range(0, igrid):
            for j in range(0, jgrid):
                xr = r[j] * np.sin(theta[i])
                xx[i, j] = xr
                yr = r[j] * np.cos(theta[i])
                yy[i, j] = yr

        x = xx[:, :]
        y = yy[:, :]

        #print "xx.shape", xx.shape
        #print "Ma.shape", Ma.shape
        #print "ene.shape:", ene.shape

        '''Intend to plot v1 and density at different log(loc) ... '''

        '''Find density'''
        den = den/(2*10e27) # solar mass per unit volume to g/cm^3
        X, Y = np.meshgrid(x, y)
        Z = np.log10(den)
        ratio_1 = 3.0856e16 # kpc to km
        ratio_2 = 3.15e16 # Gyr to s
        Z1 = v1 * ratio_1 / ratio_2

        '''Density plot'''
        im = axs[num_of_pics].pcolormesh(x, y, Z, cmap = 'rainbow')
        axs[num_of_pics].set_title(filename)

        axs[num_of_pics].set_xlabel('x')
        #axs[num_of_pics].set_ylabel('y')
        #plt.title('$\\rho$ (g cm⁻³)')

        '''V_r plot'''
        '''
        im1 = axs[num_of_pics].pcolormesh(x, y, Z1, cmap='OrRd', vmin=-4000, vmax=4000)
        
        axs[num_of_pics].set_xlabel('x')
        
        '''

        #fig.savefig("./%s.png" % ('vr-o' + filename.split('.')[-1]))  # save fig as "Mach.pdf"

        fo.end()

axs[0].set_ylabel('y')

#axs[0].set_ylim(bottom = 0.)

fig.colorbar(im, ax=axs.tolist(), label='$\\rho$ /g cm⁻³')
#fig.colorbar(im1, ax=ax2).set_label(label='$V_{r}/km s⁻¹$')
plt.show()