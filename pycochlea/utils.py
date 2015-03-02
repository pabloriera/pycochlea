from __future__ import division
import numpy as np

ascale = 0.001

def ftoerb(f):
	return 24.7 * (4.37 * f/1000 + 1)

def ftoerbscale(f):
    return 21.4*np.log10(4.37*f/1000+1)

def ftocb(f):
    return 25+75*(1+1.4*(f/1000)**2)**0.69

def rms2db(Arms):
    return 20*np.log10(Arms/20e-6/ascale)

def db2rms(Idb):
    return 20e-6*10**(Idb/20)*ascale

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


def fluctuationcv(Y,fs,fband=100,olap=0.25):

    from scipy.signal import hilbert,hanning
    

    N = Y.shape[0]
    nchan = Y.shape[1]

    winsize = int(fs*1/fband)
    olap = int(winsize*olap)

    env = abs(hilbert(Y)).T
    
    senv = env.mean(0)

    if N < winsize:
        chunks = [range(N)]
    else:
        chunks = [ range(i,i+winsize) for i in range(0, N-winsize, olap) ]

    print len(chunks)

    ST0 = np.zeros(len(chunks))
    ME0 = np.zeros(len(chunks))
    ST = np.zeros((nchan,len(chunks)))
    ME = np.zeros((nchan,len(chunks)))
    
    for nc,c in enumerate(chunks):
        win = hanning(len(c))
        win /= win.sum()

        x = senv[c]
        ME0[nc] = (win*x).sum()
        ST0[nc] = (win*(x-ME0[nc])**2).sum()
        

        x = env[:,c] 

        ME[:,nc] = (win*x).sum(1)
        ST[:,nc] = (win*( (x.T - ME[:,nc]).T )**2).sum(1)

    S0 = np.mean(ST0)
    M0 = np.mean(ME0)
    
    S = np.mean(ST,axis=1)
    M = np.mean(ME,axis=1)
    
    R = np.mean(ST/ME,axis=1)

    return S0,M0,S,M,R


def fluctuationcv_(env,fs,fband=100,olap=0.25):

    from scipy.signal import hilbert

    N = env.size

    winsize = int(fs*1/fband)
    olap = int(winsize*olap)
    
    if N < winsize:
        chunks = [range(N)]
    else:
        chunks = [ range(i,i+winsize) for i in range(0, N-winsize, olap) ]

    ST = np.zeros(len(chunks))
    ME = np.zeros(len(chunks))
    
    for nc,c in enumerate(chunks):
        win = hamming(len(c))
        win /= win.sum()
        x = env[c]
        ME[nc] = (win*x).sum()
        ST[nc] = (win*(x-x.mean())**2).sum()
        

    ST = np.sqrt(ST)
    S = np.mean(ST)
    M = np.mean(ME)
    
    R = np.mean(ST/ME)

    return S,M,R

def modulation_filter_bank(Y,fs,nf=10,fmin=2,fmax=100):
    
    from scipy.signal import lfilter
    
    fcs = np.logspace(np.log10(fmin), np.log10(fmax), nf )

    X = np.zeros((Y.shape[0],fcs.size))
    for i,fc in enumerate(fcs):
        if fc>10:
            B = fc/2
        else:
            B = 5
        
        R = np.exp(-np.pi*B/fs)
        c = np.exp(1j*2*np.pi*fc/fs)
        b = [(1-R)]
        a = [1,-c*R]
        
        X[:,i] = np.sqrt(  (lfilter(b,a,Y).real **2 ).mean(1) )
        
    return fcs,X

def modulation_filter_bank_full(Y,fs,nf=10,fmin=2,fmax=100,gpu=False):
    
    from scipy.signal import lfilter
    from pycochlea import gpuFilterbank
    
    fcs = np.logspace(np.log10(fmin), np.log10(fmax), nf )

    X = {}
    for i,fc in enumerate(fcs):
        if fc>10:
            B = fc/2
        else:
            B = 5
        
        R = np.exp(-np.pi*B/fs)
        c = np.exp(1j*2*np.pi*fc/fs)
        b = [(1-R)]
        a = [1,-c*R]
        
        if gpu: #not working with complex filter   
            gpuFilter = gpuFilterbank(newb,newa, nchan)

            gpuFilter.set_input(Y)
            X[i] = gpuFilter.process()

        else:
            X[i]= lfilter(b,a,Y).real

        return fcs,X