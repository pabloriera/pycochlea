from __future__ import division
import numpy as np
from cochlea_model import *

def ftoerb(f):
	return 24.7 * (4.37 * f/1000 + 1)

def ftoerbscale(f):
    return 21.4*np.log10(4.37*f/1000+1)


def db2rms(Idb):
    return 20e-6*10**(Idb/20)*0.0001

def normalize2db(signal,Idb):
    rms = np.sqrt(np.mean(signal**2));
    signal = db2rms(Idb)*signal/rms;
    return signal

def cosineramp(dur,N,fs=44100):

    ramp = np.ones(N)
    t = np.arange(N)/fs
    ramp[t<=dur] = 0.5 * (1 + np.cos(np.pi/dur * (t[t<=dur] - dur) ))
    return ramp


def tukeywin(window_length, alpha=0.5):
    '''The Tukey window, also known as the tapered cosine window, can be regarded as a cosine lobe of width \alpha * N / 2
    that is convolved with a rectangle window of width (1 - \alpha / 2). At \alpha = 1 it becomes rectangular, and
    at \alpha = 0 it becomes a Hann window.
 
    We use the same reference as MATLAB to provide the same results in case users compare a MATLAB output to this function
    output
 
    Reference
    ---------
 
	http://www.mathworks.com/access/helpdesk/help/toolbox/signal/tukeywin.html
 
    '''
    # Special cases
    if alpha <= 0:
        return np.ones(window_length) #rectangular window
    elif alpha >= 1:
        return np.hanning(window_length)
 
    # Normal case
    x = np.linspace(0, 1, window_length)
    w = np.ones(x.shape)
 
    # first condition 0 <= x < alpha/2
    first_condition = x<alpha/2
    w[first_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[first_condition] - alpha/2) ))
 
    # second condition already taken care of
 
    # third condition 1 - alpha / 2 <= x <= 1
    third_condition = x>=(1 - alpha/2)
    w[third_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[third_condition] - 1 + alpha/2))) 
 
    return w
