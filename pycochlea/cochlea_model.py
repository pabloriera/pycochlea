from __future__ import division

import timeit
import numpy as np
from utils import *

import ctypes
from numpy.ctypeslib import ndpointer

from sys import platform as _platform

if _platform == "linux" or _platform == "linux2":
    lib = "/pycochlea.so"
elif _platform == "darwin":
    lib = "/pycochlea.so"
elif _platform == "win32":
    lib = "/pycochlea.dll"

cochlea = ctypes.cdll.LoadLibrary(__file__[:-17]+lib).cochlea

cochlea.restype = None

#cochlea(flotante *X_t, flotante *stimulus, flotante* weight, flotante* sdivm, flotante* ddivm, flotante* params, int* dimension, flotante* solver_opts ) 
#params dv, dvy2, alpha, beta
#dimensions n_t, n_osc,fs
#solver_opts solver, abs_tol, rel_tol 

float_t = ctypes.c_double

cochlea.argtypes = [ndpointer(float_t , flags="C_CONTIGUOUS"),
                    ndpointer(float_t , flags="C_CONTIGUOUS"),
                    ndpointer(float_t , flags="C_CONTIGUOUS"),
                    ndpointer(float_t , flags="C_CONTIGUOUS"),
                    ndpointer(float_t , flags="C_CONTIGUOUS"),
                    ndpointer(float_t , flags="C_CONTIGUOUS"),
                    ndpointer(float_t , flags="C_CONTIGUOUS"),
                    ndpointer(ctypes.c_int , flags="C_CONTIGUOUS"),
                    ndpointer(float_t , flags="C_CONTIGUOUS")]

def _pre_cochlea( signal, weight, ww , dd , params_dic, dimensions_dic ,solver_dic):
    
    flotante = np.float64
    
    dimensions_dic['n_t'] = len(signal)

    params_order = ['alpha','beta','gamma','delta','epsilon','zeta','eta','nu','fluid']
    dimensions_order = ['n_t','n_osc','fs','dec']
    solver_opts_order = ['solver','abs_tol','rel_tol']

    params = np.array([params_dic[ll] for ll in params_order]).astype(flotante)
    dimensions = np.array([dimensions_dic[ll] for ll in dimensions_order]).astype(np.int32)
    solver_opts = np.array([solver_dic[ll] for ll in solver_opts_order]).astype(flotante)

    n_t = dimensions_dic['n_t']
    n_osc = dimensions_dic['n_osc']
    dec = dimensions_dic['dec']
    n_t_dec = int(n_t/dec)+n_t%dec


    tt = np.zeros( (n_t_dec, 1) ).astype(flotante)
    X_t = np.zeros( (n_t_dec, n_osc*2) ).astype(flotante)
    signal = signal.astype(flotante)
    weight = weight.astype(flotante)
    ww = ww.astype(flotante)
    dd = dd.astype(flotante)
    
    return X_t,tt,  signal, weight, ww , dd , params, dimensions, solver_opts
    
def _cochlea( signal, weight, ww , dd , params_dic, dimensions_dic ,solver_dic):

    n_t = dimensions_dic['n_t']
    n_osc = dimensions_dic['n_osc']
    fs = dimensions_dic['fs']
    
    X_t,tt,  signal, weight, ww , dd , params, dimensions, solver_opts = _pre_cochlea( signal, weight, ww , dd , params_dic, dimensions_dic ,solver_dic)

    tic = timeit.default_timer()
    cochlea(X_t,tt,  signal, weight, ww , dd , params, dimensions, solver_opts)
    toc = timeit.default_timer()
    
    
    Y = X_t[:,:n_osc]
    V = X_t[:,n_osc-1:-1]
    
    if np.sum(np.isnan(Y))!=0:
        print 'Simulation exploted'        
    
    return {'Y':Y,'V':V,'tt':tt,'tictoc':(toc-tic)/n_t*fs}

def time_cochlea( signal, weight, ww , dd , params_dic, dimensions_dic ,solver_dic,nrep ):
    
    n_t = dimensions_dic['n_t']
    fs = dimensions_dic['fs']
    
    X_t,tt,  signal, weight, ww , dd , params, dimensions, solver_opts = _pre_cochlea( signal, weight, ww , dd , params_dic, dimensions_dic ,solver_dic)
    
    t = 0
    for i in xrange(nrep):
        print i
        tic = timeit.default_timer()
        cochlea(X_t,tt,  signal, weight, ww , dd , params, dimensions, solver_opts)
        toc = timeit.default_timer()
        t += (toc-tic)/n_t*fs
    
    return t/nrep
    
        
def phase_response(Y,method='bins'):
    from pylab import find

    n_t,nchan = Y.shape
    
    if method=='fft':
       
        FY = np.fft.fft(Y.T)   
        
        FY_db = 20*np.log10(np.abs(FY))
        chf = np.argmax(np.max(FY_db,axis=0))
        #ph = unwrap(angle(FY),discont=0.1)
        ph = np.angle(FY)
        x = ph[:,chf].copy()
        
        k=1
        for i in xrange(1,nchan-1):
            if x[i]>x[i-1]:
                x[i]-=k*2*np.pi
                if (x[i]-x[i-1])>0.2:
                    k+=1
                    x[i]-=2*np.pi
        
        return (x-x[0])/np.pi/2
        
    elif method=='bins':
        
        DY = np.diff(Y,axis=0)
        F = np.logical_and(DY[:-1]>0,  DY[1:]<0)
        F = F[n_t*1/2:,:]
        xp = np.zeros(nchan)
        xp[0]=find(F[:,0])[0]
        
        for i in xrange(1,nchan):
            if find(F[xp[i-1]:,i]).size>0:
                xp[i]=xp[i-1]+find(F[xp[i-1]:,i])[0]
            
        return -(xp-xp[0])
        
def middle_ear(S,fs):
    
    from scipy.signal import zpk2tf,bilinear,freqs,lfilter
    
    # p1 = -800 + 2100j
    # p0 = -250 + 680j
    # p2 = -2000 + 10000j
    # z0 = -1800
    # z1 =  (-000 + 5000j)*2*np.pi
    # z2 = -1000000 + 20000*2*np.pi*1j
    # zeros = np.array([z0, z0,z1,conj(z1),z2,z2 ])
    # poles = np.array([p0, conj(p0), p1, conj(p1),p2,conj(p2)] )

    z0 = -800
    p0 = -250 + 680j
    p1 = -3000 + 4000j
    zeros = np.array([z0, z0])
    poles = np.array([p0, np.conj(p0), p1, np.conj(p1)] )
    # use an arbitrary gain here, will be normalized afterwards
    b, a = zpk2tf(zeros, poles * 2 * np.pi, 1.5e9)
    # normalize the response at 1000Hz (of the analog filter)
    resp = np.abs(freqs(b, a, [1000*2*np.pi])[1])  # response magnitude
    b /= resp
    bd, ad = bilinear(b, a, fs)
    
    return lfilter(bd,ad,S)

def runcochlea(signal,data):

    dimensions={}
    params={}
    solver={}

    if 'mef' not in data.keys():
        data['mef']=10
    
    dimensions['n_osc'] = data['nchan']
    dimensions['n_t'] = len(signal)
    dimensions['fs'] = data['fs']
    
    dx = data['length']/data['nchan']
       
    kappa = data['mass']*data['height']/2/data['density']
    
    params['alpha'] = kappa / dx / data['height']*data['mef']
    params['beta'] = dx / data['height']*data['mef']
    params['gamma'] = data['gamma']
    params['delta'] = data['delta']
    params['epsilon'] = data['alphita']
    params['zeta'] = -data['gammita']*data['alphita']
    params['eta'] = data['eta']
    params['nu'] = data['nu']
    params['fluid'] = data['fluid']
    
    if 'dec' not in data.keys():
        dimensions['dec']=1
    else:
        dimensions['dec']=data['dec']
    
    if 'solver' in data.keys():
        solver['solver'] = data['solver']
        solver['abs_tol'] = data['abs_tol']
        solver['rel_tol'] = data['rel_tol']
    
    else:
        solver['solver'] = 0
        solver['abs_tol'] = 1e-5
        solver['rel_tol'] = 1e-5
    
    weight = np.zeros(data['nchan'])
    weight[0] = 1
    
    ff = np.flipud( np.logspace(np.log10(data['fmin']),np.log10(data['fmax']),data['nchan']))
    
    w = 2*np.pi*ff
    ww = w**2
    dd = w/data['Q']    
    
    if data['middle_ear']:
        signal = middle_ear(signal,data['fs'])
        
    #signal = data['mass']*np.convolve(signal,[1,-2,1],mode='same')*data['fs']**2
    signal = data['mass']*np.convolve(signal,[1,-2,1],mode='same')*data['fs']**2
    
    X = _cochlea( signal, weight, ww , dd , params, dimensions ,solver )
    
    return X



    
def numpy_cochlea(pure_tone,data):
    
    def blockset(A,n,b):
        N = A.shape[0]
        a = n % 2 + (n//2)*N
        A.flat[a::N*2+2]=b
        return A
    
    #Runge-Kutta 4 Integrator
    def rk4(t, dt, x, f):
        k1 = dt * f(t, x)
        k2 = dt * f(t + 0.5*dt, x + 0.5*k1)
        k3 = dt * f(t + 0.5*dt, x + 0.5*k2)
        k4 = dt * f(t + dt, x + k3)
        return t + dt, x + (k1 + 2*(k2 + k3) + k4)/6.0     
    
    
    def f1(x): return x
    
    N = data['nchan']
    fs = data['fs']
    
    dx = data['length']/data['nchan']
       
    kappa = data['mass']*data['height']/2/data['density']
    
    alpha = kappa / dx / data['height']  
    beta = dx / data['height']    
    gamma = data['gamma']
    delta = data['delta']
    epsilon = data['alphita']
    zeta = data['alphita']*data['gammita']
    eta = data['eta']
    fluid = data['fluid']
    
    weight = np.zeros(data['nchan'])
    weight[0] = 1
    
    ff = np.flipud( np.logspace(np.log10(data['fmin']),np.log10(data['fmax']),data['nchan']))
    
    w = 2*np.pi*ff
    ww = -w**2
    dd = -w/data['Q']   
    
    #Prepare fluid matrix coefs
    a = alpha/beta * np.ones(N-1)
    b = -2*alpha/beta-1 * np.ones(N)
    c = alpha/beta * np.ones(N-1)
    
    b[0] = -alpha-1;
    c[0] = alpha;
     
    b[N-1] = 1.0;
    a[N-2] = 0.0;
    
    N_delta = N*delta
    k_1  = int(np.floor(N_delta ))
    k_2  = k_1 + 1
    ee = N_delta - k_1
    
    #lateral feed matrix
    IL = -zeta*( np.diag(np.ones(N-k_1),-k_1)*(1-ee) + np.diag(np.ones(N-k_2),-k_2)*ee) \
                + epsilon*(np.diag(np.ones(N-k_1),k_1)*(1-ee) + np.diag(np.ones(N-k_2),k_2)*ee)
    
    #Fluid matrix
    P = np.diag(a,-1)+np.diag(b)+np.diag(c,1)
    IP = np.linalg.inv(P)
#    IP[:,-1]=0
    
    G = np.zeros((N*2,N*2))
    C = np.zeros((N*2))
    D = np.zeros((N*2,N*2))
    D2 = np.zeros((N*2,N*2))
    L = np.zeros((N*2,N*2))
    G = blockset(G,1, np.ones(N))
    G = blockset(G,2, ww)
    G = blockset(G,3, dd)
    
    #Input stimulus matrix weights
    C[1] = 1;
    
    for i in xrange(N):
        D[i*2+1,:-2:2]  = IP[i,:-1]*ww[:-1]
        D[i*2+1,1:-1:2] = IP[i,:-1]*dd[:-1]
    
    for i in xrange(N):
#    D2[i*2+1,::2]  = IP[i,:]
        D2[i*2+1,1::2] = IP[i,:]
    
    for i in xrange(k_1,N-k_1):
        L[i*2,:-2:2]  = IL[i,:-1]        
        
    #main matrix
    A = (D+G).dot(np.eye(2*N)+L)
        
    #The system linear coefficients
    
    #Creates B matrix, for the nonlinear functioni f1
    #B = np.zeros((N*2,N*2))
    #B = blockset(B,3, np.ones(N))
    #B = np.float32(B)
    
    A_db = pure_tone['amplitude_db']
    Amp = 20e-6*10**(A_db/20)
    f0 = pure_tone['f0']
    w0 = 2*np.pi*f0
     
    def env(t,T,r):
        if t<r/2:
            return 0.5 * (1 + np.cos(2*np.pi/r * (t - r/2) ))
        elif t>T-r/2:
            return 0.5 * (1 + np.cos(2*np.pi/r * (t - T + r/2) ))
        else:
            return 1
    
    dt = 1/fs
    t = 0
    dur = pure_tone['duration']
    n_t = int(dur/dt)
    
    def signal0(t): return Amp*w0**2*np.sin(w0*t)*env(t,dur,0.2)
        
    def dfun(t,x): return A.dot(x) + D2.dot(C*signal0(t))# + B.dot(f1(x))
    
    x = np.zeros(N*2)
    out = np.zeros((N*2,n_t))        
    
    ii = -1    
    for i in xrange(n_t):
        ii+=1
        t, x = rk4(t,dt,x,dfun)            
        out[:,i] = x

    X = {'Y':out[::2,:],'V':out[1::2,:]}
    
    return X

def pures2cochlea(pure_tone,data):       
     
    n_t = np.ceil(data['fs']*pure_tone['duration'])
    t = np.arange(n_t)/data['fs']
        
    signal = np.zeros(t.shape)  
    
    for a,f in zip(pure_tone['amplitude_db'],pure_tone['f0']):
        
        A = db2rms(a)
        signal += A*np.sin(2*np.pi*f*t)    
    
    signal = tukeywin(n_t,0.2)*signal
     
    X = runcochlea(signal,data)
    
    return X
    
def pure2cochlea(pure_tone,data):       
    
    n_t = np.ceil(data['fs']*pure_tone['duration'])
    t = np.arange(n_t)/data['fs']

    A = db2rms(pure_tone['amplitude_db'])
    signal = tukeywin(n_t,0.2)*A*np.sin(2*np.pi*pure_tone['f0']*t)    

    X = runcochlea(signal,data)
    
    return X

