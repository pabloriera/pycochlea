# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 12:45:40 2015

@author: miles
"""

import numpy as np
from utils import *
from cochlea_model import *

from brian.hears import LinearFilterbank
from brian.stdunits import kHz, Hz, ms

def LowPass_filter_(fc,fs):
    
    TWOPI = 2*np.pi
    
    c = 2.0 * fs
    c1LP = ( c - TWOPI*fc ) / ( c + TWOPI*fc )
    c2LP = TWOPI*fc / (TWOPI*fc + c)
    
    b = np.array([c2LP,c2LP])
    a = np.array([1,-c1LP])
    
    return b,a
    
    
def filter_coeff_one2n(a,n):
    import sympy
    import numpy as np
    
    x = sympy.symbols("x")
    
    formula = (a[0] + a[1]*x) ** n
    newa = formula.expand().as_poly().all_coeffs()
        
    return np.array(np.float64(newa[::-1]))


class LowPass_filter(LinearFilterbank):

    def __init__(self,source,cf,fc,gain,order):
        nch = len(cf)
        TWOPI = 2*np.pi
        self.samplerate =  source.samplerate
        c = 2.0 * self.samplerate
        c1LP = ( c/Hz - TWOPI*fc ) / ( c/Hz + TWOPI*fc )
        c2LP = TWOPI*fc/Hz / (TWOPI*fc + c/Hz)
        
        b_temp = np.array([c2LP,c2LP])
        a_temp = np.array([1,-c1LP])
        
        filt_b = np.tile(b_temp.reshape([2,1]),[nch,1,order])               
        filt_a = np.tile(a_temp.reshape([2,1]),[nch,1,order]) 
        filt_b[:,:,0] = filt_b[:,:,0]*gain

        LinearFilterbank.__init__(self, source, filt_b, filt_a)
        
def IHC_transduction(x,slope,asym,sign): 
    corner    = 80
    strength  = 20.0e6/10**(corner/20)     
    xx = sign*np.log(1.0+strength*abs(x))*slope
    ind = x<0
    splx   = 20*np.log10(-x[ind]/20e-6);
    asym_t = asym-(asym-1)/(1+np.exp(splx/5.0));
    xx[ind] = -1/asym_t*xx[ind]
    return xx

def zilany2014_run_synapse(vihc,ff,args,channels):
    from cochlea.zilany2014 import _zilany2014

    fs = args['fs']
    
    powerlaw = args['powerlaw']
    anf_num = args['anf_num']
    ffGn = args['ffGn']

    duration = vihc.shape[0] / fs
    anf_types = np.repeat(['hsr', 'msr', 'lsr'], anf_num)

    nested = []
    for i in channels:

        cf = ff[i]
        vihc_canal = vihc[:,i]

        synout = {}
        trains = []
        for anf_type in anf_types:
    
            if anf_type not in synout:
                ### Run synapse
                synout[anf_type] = _zilany2014.run_synapse(
                    fs=fs,
                    vihc=vihc_canal,
                    cf=cf,
                    anf_type=anf_type,
                    powerlaw=powerlaw,
                    ffGn=ffGn
                )
    
            ### Run spike generator
            spikes = _zilany2014.run_spike_generator(
                synout=synout[anf_type],
                fs=fs,
            )
    
    
            trains.append({
                'spikes': spikes,
                'duration': duration,
                'cf': cf,
                'type': anf_type
            })
            
        nested.append(trains)
            
    return nested

def synapse(vihc,ff,args,channels):
    from cochlea.zilany2014 import _zilany2014

    fs = args['fs']
    
    powerlaw = args['powerlaw']
    ffGn = args['ffGn']

    duration = vihc.shape[0] / fs
    anf_types = ['hsr', 'msr', 'lsr']

    synout = {}
    
    for anf_type in anf_types:
        synout[anf_type ] = np.zeros((len(channels),vihc.shape[0]))   

    for i in channels:

        cf = ff[i]
        vihc_canal = vihc[:,i]
       
        for anf_type in anf_types:
    
            syn = _zilany2014.run_synapse(
                fs=fs,
                vihc=vihc_canal,
                cf=cf,
                anf_type=anf_type,
                powerlaw=powerlaw,
                ffGn=ffGn
            )
        
            synout[anf_type ][i,:] = syn
            
    return synout

def ihc(Y,fs,Yscale=1e2,fcut=3000):
    from scipy.signal import lfilter

    blop,alop = LowPass_filter_(fcut,fs)
    alop = filter_coeff_one2n(alop,7)
    blop = filter_coeff_one2n(blop,7)

    Yihc = IHC_transduction(Y*Yscale,slope = 0.1,asym = 3.0,sign=1)
    vihc = lfilter(blop,alop,Yihc.T).T
    return vihc

def ihc_synapse(Y,data,Yscale=0.7e2,fcut=3000,vihc_scale=1.5,channels=0):

    fs = int( data['fs']/data['dec'] )

    vihc = ihc(Y,fs,Yscale=Yscale,fcut=fcut)*vihc_scale

    ff = np.flipud( np.logspace(np.log10(data['fmin']),np.log10(data['fmax']),data['nchan']))

    args = { 'fs': fs,  'powerlaw': 'approximate','ffGn': False}

    if channels==0:
        channels = xrange(data['nchan'])

    synout = synapse(vihc,ff,args,channels)

    return synout,vihc

def ihcan_synapse(Y,data,Yscale=1e2,vihc_scale=2,anf_num = (100,50,30),channels=0):
    
    import itertools
    import pandas as pd
    from brian.hears import Sound, FunctionFilterbank

    ff = np.flipud( np.logspace(np.log10(data['fmin']),np.log10(data['fmax']),data['nchan']))

    pre_ihc = Sound(Y*Yscale,samplerate=data['fs']*Hz)

    pre_ihc2 = FunctionFilterbank(pre_ihc,IHC_transduction,slope = 0.1,asym = 3.0,sign=1)

    vihc = LowPass_filter(pre_ihc2,ff,3000,1.0,7)*vihc_scale

    args = { 'fs': data['fs'], 'anf_num': anf_num,  'powerlaw': 'approximate','ffGn': False}

    if channels==0:
        channels = xrange(data['nchan'])
    
    nested = zilany2014_run_synapse(vihc.process(),ff,args,channels)
    
    trains = itertools.chain(*nested)
    anf = pd.DataFrame(list(trains))

    return anf

def pure2periphery(pure_tone,data,channels=0):       
    
    X = pure2cochlea(pure_tone,data)
    Y = X['Y']
    syn,vihc = ihc_synapse(Y,data,channels=channels)

    return {'Y':Y,'syn':syn}

def pure2periphery_(pure_tone,data,periphery_data={}):       
    
    X = pure2cochlea(pure_tone,data)
    Y = X['Y']
    anf_num = (50,30,20)
    channels = periphery_data['channels']
    anf = ihcan_synapse(Y,data,Yscale=1e2,vihc_scale=2,anf_num =anf_num,channels=channels)

    return {'Y':Y,'anf':anf}
