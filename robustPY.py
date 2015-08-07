# This demonstrates a basic implementation of the
# ROBUST_f2py.F90 wrapped code that imports 'robust.so'
#
# Note: To get around allocatable arrays in f2py, this version returns 'nsteps'
# which is the number of tc,tv incrememnts computed inside the module.
# you can use this to trim the returned arrays from their default size as shown.

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mtb
import pylab as py
import sys
import robust as tctv

plt.ion()
plt.close()


params = {'legend.fontsize': 19}
py.rcParams.update(params)

def drawfig(mstar,tc,tv):
    fig, axs = plt.subplots(1,1,figsize=(8, 5), facecolor='w', edgecolor='k')
    plt.subplots_adjust(hspace=0.2, bottom=0.1, wspace=0,left=0.09,right=0.98,top=0.99)

    axs.plot(mstar,tc,'-',color='r',zorder=2,linewidth=1,label=r'$\mathsf{T_c}$')
    axs.scatter(mstar, tc,color='r')

    axs.plot(mstar,tv,'-',color='k',zorder=2,linewidth=1,label=r'$\mathsf{T_v}$')
    axs.scatter(mstar, tv,color='k')

    axs.set_xlim([13.5,21.5])
    axs.set_ylim([-6,6])
    axs.axhline(y=-3, xmin=0, xmax=1,linewidth=1, color = 'g',linestyle='--')
    axs.axhline(y=3, xmin=0, xmax=1,linewidth=1, color = 'g',linestyle='--')
    axs.axvline(x=20, ymin=-3, ymax=3,linewidth=1, color = 'b',linestyle='--')
    axs.set_xlabel(r'$\mathsf{Trial\ apparent\ magnitude\ (m_*)}$',fontsize=16)
    axs.set_ylabel(r'$\mathsf{Completness\ limits}$',fontsize=16)

    lg = axs.legend(loc='upper left',ncol=1,fontsize=14)
    lg.draw_frame(False)
    
    plt.show()

    return


print tctv.robust.tctv_r01.__doc__
print


## READ IN THE TEST DATA ##
mt,am,mu = np.loadtxt('testdata.txt',usecols=(0,1,2), unpack='true')


# SET INPUT PARAMETERS.
ngal=len(mt)
mag_min=np.min(mt)
mag_max=np.max(mt)
bin=0.1


###########  - MAIN - #################################
#calling the R01 version of ROBUST:

tc,tv,mstar,nsteps=tctv.robust.tctv_r01(mag_min,mag_max,mu,mt,am,bin,ngal)

tc= tc[0:nsteps]
tv= tv[0:nsteps]
mstar= mstar[0:nsteps]

#plot the output
drawfig(mstar,tc,tv)


#calling the JTH07 version of ROBUST:
#setting the delta_mu,am width here.
delta_am=0.5
delta_mu=0.5
tc,tv,mstar,nsteps=tctv.robust.tctv_jth07(mag_min,mag_max,mu,mt,am,bin,delta_am,delta_mu,ngal)


tc= tc[0:nsteps]
tv= tv[0:nsteps]
mstar= mstar[0:nsteps]

drawfig(mstar,tc,tv)
