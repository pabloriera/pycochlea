# -*- coding: utf-8 -*-
from __future__ import division
import pylab as pl
import numpy as np
from pycochlea import pure2cochlea_explore_multiprocess


#%%
  
f0 = 1000
A_db_list = np.arange(1,11)*10
f0_list = [f0]*len(A_db_list)
duration_list = [0.2]*len(A_db_list)

pure_tone = {'amplitude_db':A_db_list,'f0':f0_list,'duration':duration_list}


data = {'fs': 100000.0,'fmin': 20.0, 'fmax': 20000, 'height': 0.001, 'Q': 5, 
        'delta': 0.004, 'middle_ear': True, 'nchan': 500, 
        'density': 1000, 'length': 0.035,
        'mass': 0.03, 'gamma': 2, 'nu':0, 'gammita': 0.3,'alphita': 0.15,
        'eta':1*1e11,'fluid':1,'solver':0,'abs_tol':1e-8,'rel_tol': 1e-8}


re = pure2cochlea_explore_multiprocess(pure_tone,data,nprocs = 4)
#%%

ff = np.flipud( np.logspace(np.log10(data['fmin']),np.log10(data['fmax']),data['nchan']))

cmap = pl.cm.get_cmap(name='spectral')
cmap = [cmap(i) for i in np.linspace(0, 0.9, len(re))]

A_db_list = [r['amplitude_db'] for r in re]
A_db_list.sort()

V_db_max = np.zeros(len(A_db_list))
Y_max = np.zeros(len(A_db_list))
V_db = np.zeros((len(A_db_list),data['nchan']))
Y_db = np.zeros((len(A_db_list),data['nchan']))
V_mean = np.zeros((len(A_db_list),data['nchan']))
freqs = [f0*0.95,f0,f0*1.2,f0*4]

for i,r in enumerate(re):
    
    j = A_db_list.index(r['amplitude_db'])
    Y = r['out']  
    Y_db[j,:] = 20*np.log10(np.sqrt(np.mean(Y[Y.shape[0]*3/4:,:]**2,axis=0))/1e-6)
    

#%%
num=12198
pl.figure(num=num,figsize=(8,9))
pl.clf()

for j,A in enumerate(A_db_list):
    print A,Y_max[j]
    pl.subplot(2,1,1)
    pl.semilogx(ff,Y_db[j,:],'-',label=A,color=cmap[j])
    pl.subplot(2,1,2)
    pl.semilogx(ff,Y_db[j,:]-Y_db[j,0],'-',label=A,color=cmap[j])

pl.subplot(2,1,1)
pl.ylabel('Desplazamientos RMS MB \n (dB re $1\mu m$)')
pl.xlabel(u'Frecuencia característica (Hz)')
pl.xlim(ff[0],ff[-1])
pl.ylim(-90,35)

for f in freqs:
    pl.gca().arrow(f, 35, 0, -5, head_width=1, head_length=2, fc='k', ec='k')

pl.subplot(2,1,2)
pl.ylabel('Ganancia MB (dB)')
pl.xlabel(u'Frecuencia característica (Hz)')
pl.xlim(ff[0],ff[-1])
pl.ylim(0,65)
pl.legend(prop={'size':15},title='$dB_{SPL}$')

pl.tight_layout()

pl.figure(12199)
pl.cla()
A_db_list = np.array(A_db_list)

ixs = []

for f in freqs:
    ixs.append(np.argwhere(f>ff)[0])

Y_aux = Y_db[:,ixs]

pl.plot(A_db_list,A_db_list-120,'k--',lw=3,alpha=0.5)

for j in range(Y_aux.shape[1]):
    pl.plot(A_db_list,Y_aux[:,j],'.-',lw=2,label=np.around(freqs[j]) )

pl.legend(loc=4,prop={'size':15},title='$F_c (Hz)$')
pl.ylabel('Desplazamientos RMS MB \n (dB re $1\mu m$)')
pl.xlabel(u'Intensidad estímulo ($dB_{SPL}$)')
pl.ylim(-150,40)
pl.tight_layout()