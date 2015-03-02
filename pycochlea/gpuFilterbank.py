# -*- coding: utf-8 -*-
from __future__ import division
"""
Created on Wed Aug 20 13:39:22 2014

@author: miles
"""

#%%

import numpy as np
import pycuda.autoinit
#import pycuda.driver as drv
from pycuda.compiler import SourceModule
from pycuda import gpuarray
import re

'''
``precision``
    Should be single for older GPUs.
``unroll_filterorder``
    Whether or not to unroll the loop for the filter order, normally
    this should be done for small filter orders. By default, it is done
    if the filter order is less than or equal to 32.
    
    
    The filterbank is implemented as a direct II transposed structure.
This means that for a single channel and element of the filter cascade,
the output y for an signal x is defined by::

    a[0]*y[m] = b[0]*x[m] + b[1]*x[m-1] + ... + b[m]*x[0]
                          - a[1]*y[m-1] - ... - a[m]*y[0]

using the following difference equations::

    y[i] = b[0]*x[i] + z[0,i-1]
    z[0,i] = b[1]*x[i] + z[1,i-1] - a[1]*y[i]
    ...
    z[m-3,i] = b[m-2]*x[i] + z[m-2,i-1] - a[m-2]*y[i]
    z[m-2,i] = b[m-1]*x[i] - a[m-1]*y[i]

where i is the output sample number.

    y[n] = b0 * x[n] + z[0,n]
    z[0,n] = b1 * x[n] - a1 * y[n] + z[1,n]
    z[1,n] = b2 * x[n] - a2 * y[n] 
    
'''

class gpuFilterbank():

    def __init__(self, b, a, nchan = 1, precision='double', unroll_filterorder=None):
               
        a = np.array(a)
        b = np.array(b)
        
        self.nchan = nchan
        
        if b.shape[0]==1 or b.ndim == 1:
            b = np.tile(b,(nchan,1))
            a = np.tile(a,(nchan,1))
            
        if precision=='double':
            self.pd = self.precision_dtype=np.float64
        else:
            self.pd = self.precision_dtype=np.float32
            
        #channels, order
        n,m = b.shape
        
        self.filt_b_gpu = np.array(b, dtype=self.precision_dtype)
        self.filt_a_gpu = np.array(a, dtype=self.precision_dtype)
        self.filt_b_gpu = gpuarray.to_gpu(self.filt_b_gpu.T.flatten()) # transform to Fortran order for better GPU mem
        self.filt_a_gpu = gpuarray.to_gpu(self.filt_a_gpu.T.flatten()) # access speeds
            
        unroll_filterorder = unroll_filterorder
        
        if unroll_filterorder is None:
            if m<=32:
                unroll_filterorder = True
            else:
                unroll_filterorder = False        
        
        code='''
        #include "stdio.h"
        #define y(s,j) y[(j)+nchan*(s)]
        #define x(s,j) x[(j)+nchan*(s)]
        #define a(i,j) a[(j)+nchan*(i)]
        #define b(i,j) b[(j)+nchan*(i)]
       
        __global__ void filt(SCALAR *b, SCALAR *a, SCALAR *x, SCALAR *y, int numsamples)
        {
            int j = blockIdx.x * blockDim.x + threadIdx.x;
            if(j>=nchan) return;
            
            SCALAR z[m];
            
            //if(j==0)printf("%g\\t%g\\n",b[0],b[nchan]);
            
            //if(j==0)printf("%g\\t%g\\n",a[0],a[nchan]);
            
            for(int i=0; i<m-1; i++)
                z[i] = 0;
                
            for(int s=0; s<numsamples; s++)
            {
                //if(j==0)printf("%g ",x[s*nchan]);
        '''
        loopcode='''
        y(s,j) = b(0,j)*x(s,j) + z[0];
        '''
        if unroll_filterorder:
            for i in range(m-2):
                loopcode+=re.sub('\\bi\\b', str(i), '''
                z[i] = b(i+1,j)*x(s,j) + z[i+1] - a(i+1,j)*y(s,j);
                ''')
        else:
            loopcode+='''
            for(int i=0;i<m-2;i++)
                z[i] = b(i+1,j)*x(s,j) + z[i+1] - a(i+1,j)*y(s,j);
            '''
        loopcode+='''
        z[m-2] = b(m-1,j)*x(s,j) - a(m-1,j)*y(s,j);
        '''
        
        code+=loopcode
        code+='''
            }
        }
        '''
        
        code=code.replace('SCALAR', precision)
        code=re.sub("\\bm\\b", str(m), code)
        code=re.sub("\\bnchan\\b", str(nchan), code)
    #    print code
        
        gpu_mod = SourceModule(code)
        self.gpu_filt_func = gpu_mod.get_function("filt")
        
        THREAD_NUM = 32
        BLOCK_NUM = int(np.ceil(nchan/THREAD_NUM))
    
        self.gridShape = (BLOCK_NUM,1,1)
        self.blockShape = (THREAD_NUM,1,1)
        
        self.gpu_filt_func.prepare((np.intp, np.intp, np.intp, np.intp, np.int32))
        
        self.setted = False
        
    def set_input(self, signal):
        
        """if signal is a numpy array its shape should be (nchans,nsamples)
           if signal is a gpuarray its shape should be (nchans*nsamples,) as
           it was created like this signal_gpu = gpuarray.to_gpu( signal.T.flatten() )
        """
        self.signal = signal
        
        if type(signal)==np.ndarray:       
            if self.nchan == signal.shape[0]:
            
                self.nsamples = np.int32(signal.shape[1])
                self.filt_x_gpu = gpuarray.to_gpu(signal.T.flatten().astype(self.pd))
                self.filt_y_gpu = gpuarray.empty_like(self.filt_x_gpu)
                self._desiredshape = signal.T.shape
                self.setted = True
                    
            else:
                print "signal must have shape as (nchan,nsamples)"
            
        elif type(signal)==gpuarray.GPUArray:                
                
                self.nsamples = np.int32( signal.size / self.nchan )
                self.filt_x_gpu = signal
                self.filt_y_gpu = gpuarray.empty_like(self.filt_x_gpu)
                self._desiredshape = (self.nsamples,self.nchan)
                self.setted = True
             
    def process(self, sync = True, ntimes = 1):
        
        if self.setted:

            self.gpu_filt_func.prepared_call(self.gridShape, self.blockShape, self.filt_b_gpu.gpudata, self.filt_a_gpu.gpudata, self.filt_x_gpu.gpudata,self.filt_y_gpu.gpudata, self.nsamples)
            
            for i in xrange(ntimes-1):

                self.filt_x_gpu.gpudata = self.filt_y_gpu.gpudata
                self.gpu_filt_func.prepared_call(self.gridShape, self.blockShape, self.filt_b_gpu.gpudata, self.filt_a_gpu.gpudata, self.filt_x_gpu.gpudata,self.filt_y_gpu.gpudata, self.nsamples)
                
            if sync:
                
                return np.reshape(self.filt_y_gpu.get(), self._desiredshape).T
                
            else:
                
                return self.filt_y_gpu
            
        else:
            print "First run set_intpu(signal)"    
