# -*- coding: utf-8 -*-
from __future__ import division
import pylab as pl
import numpy as np

from pycochlea import pures2cochlea

#%%

#pure_tones = {'amplitude_db':[A,A,A],'f0':[10000,5000,1000],'duration':0.1}
pure_tone = {'amplitude_db':[40],'f0':[1000],'duration':0.1}

data = {'fs': 100000.0, 'fmax': 20000, 'height': 0.001, 'Q': 4, 
        'delta': 0.004, 'middle_ear': False, 'gammita': 0.3,
        'nchan': 250, 'density': 1000,'nu':0, 'length': 0.035,
        'mass': 0.03, 'gamma': 2, 'fmin': 81.0, 'alphita': 0.15,
        'eta':5000*1e-8,'fluid':1}
        
data['solver'] = 0
data['abs_tol'] = 1e-8
data['rel_tol'] = 1e-8
        
ff = np.flipud( np.logspace(np.log10(data['fmin']),np.log10(data['fmax']),data['nchan']))
        
X = pures2cochlea(pure_tone,data)
print X['tictoc']

Y = X['Y']
V = X['V']

Y_db = 20*np.log10(np.max(abs(Y),axis=0))
V_db = 20*np.log10(np.sqrt(np.mean(V**2,axis=0)))
V_mean =np.sqrt(np.mean(V**2,axis=0))
Y_mean =np.sqrt(np.mean(Y**2,axis=0))

pl.figure(12354)
#pl.clf()
pl.loglog(ff,Y_mean,'-')
pl.xlim(ff[0],ff[-1])

pl.figure(12355)
pl.imshow(Y.T,aspect='auto',cmap=pl.cm.seismic)
