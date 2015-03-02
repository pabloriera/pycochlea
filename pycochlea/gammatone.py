# -*- coding: utf-8 -*-
"""
Created on Sat Nov  1 01:21:51 2014

@author: miles
"""

from brian import Hz,kHz
import numpy as np
from brian.hears import Gammatone,erbspace, Sound

def gammatone(s0,fs,cf='erb',Nf=1000):
    
    if cf=='erb':
        cf = erbspace(20*Hz, 20*kHz, Nf)
        z = 21.4*np.log10(4.37*cf/1000+1)
        return cf,z,Gammatone(Sound(s0,samplerate=fs*Hz), cf).process().T
    else:
        return Gammatone(Sound(s0,samplerate=fs*Hz), cf).process().T
    