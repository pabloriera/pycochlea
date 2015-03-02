# -*- coding: utf-8 -*-
from __future__ import division
import pylab as pl
import numpy as np
from pycochlea import synth2cochlea_explore_multiprocess


#==============================================================================
#%% 
#==============================================================================
data = {'fs': 100000.0,'fmin': 20.0, 'fmax': 20000, 'height': 0.001, 'Q': 5, 
        'delta': 0.004, 'middle_ear': True, 'nchan': 500, 
        'density': 1000, 'length': 0.035,
        'mass': 0.03, 'gamma': 2, 'nu':0, 'gammita': 0.3,'alphita': 0.15,
        'eta':1*1e11,'fluid':1,'solver':0,'abs_tol':1e-8,'rel_tol': 1e-8}

#==============================================================================
# 
#==============================================================================

def SAM(fc=1000.0,fm=70.0,m=1.0,dur=0.1,Idb=40,fs=100000):
    from pycochlea.utils import cosineramp,normalize2db
    
    t = np.arange(dur*fs)/fs
    
    signal = (1+m*np.sin(2*np.pi*fm*t))*np.sin(2*np.pi*fc*t)
    signal = signal*cosineramp(dur*0.2,t.size,fs)
    signal = normalize2db(signal,Idb)
    
    return signal
    
nsp = 10
synth_parameters = {'fm' :  np.logspace(np.log10(1), np.log10(200),nsp)}

#==============================================================================
# 
#==============================================================================

re = synth2cochlea_explore_multiprocess(SAM,synth_parameters,data,nprocs = 4)

#%%

ff = np.flipud( np.logspace(np.log10(data['fmin']),np.log10(data['fmax']),data['nchan']))

cmap = pl.cm.get_cmap(name='spectral')
cmap = [cmap(i) for i in np.linspace(0, 0.9, len(re))]

pl.figure(35654)
pl.clf()
for i,r in enumerate(re):

    Y = r['out']/1e-6 #micrometros
    j = list(synth_parameters['fm']).index(r['fm'])
    
    pl.loglog(ff,Y[Y.shape[0]*0.5:].max(axis=0),label = r['fm'],color=cmap[j] )
    pl.xlim(ff[0],ff[-1])

pl.ylim(0.001,10)