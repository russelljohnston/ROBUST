import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mtb
import pylab as py
import sys

plt.ion()
plt.close()

fin = sys.argv[1]


params = {'legend.fontsize': 19}
py.rcParams.update(params)

#fin='tctv_faint_results.out'

mstar,tc,tv = np.loadtxt(fin,usecols=(0,1,2), unpack='true')


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
raw_input()