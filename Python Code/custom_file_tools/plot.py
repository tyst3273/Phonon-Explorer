import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import h5py


#f = 'half_breathing/half_breathing_6+h-4-h_220K.hdf5'
f = 'no_couple/no_couple_6+h-4+h_220K.hdf5'
with h5py.File(f,'r') as db:
    Dim_0 = db['Dim_0'][...]
    Dim_1 = db['Dim_1'][...]
    sig = db['signal'][...]


#plt.plot(Dim_0,sig,marker='o',c='b',ms=5,lw=2)
#plt.show()

fig, ax = plt.subplots(figsize=(6,6))

im = ax.imshow(sig[...].T,cmap='viridis',origin='lower',aspect='auto',
        vmin=0,vmax=0.00016,
        #norm=LogNorm(vmin=1e-3,vmax=0.02),
        extent=[Dim_0.min(),Dim_0.max(),Dim_1.min(),Dim_1.max()],interpolation='none')
fig.colorbar(im,ax=ax,extend='both')

#ax.axis([3,9,50,100])
#ax.set_ylim([50,100])

plt.show()



