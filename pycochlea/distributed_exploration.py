from __future__ import division

from cochlea_model import *
from periphery import *
from utils import *

def chunks(L, n):
    return [L[i::n] for i in range(n)]
        
def pure2cochlea_explore_thread(pure_tone,data,nthreads=1,func=None):
    

    from threading import Thread
       
    f0_list = pure_tone['f0']
    A_db_list = pure_tone['amplitude_db']
    duration_list =pure_tone['duration']
    
    responses= []
    
    class worker(Thread):
        def __init__ (self,jj):
            Thread.__init__(self)
            self.jj = jj
            
        def run(self):
            jj = self.jj
            
            for j in jj:
                pure_tone = {'amplitude_db':A_db_list[j],'f0':f0_list[j],'duration':duration_list[j]}
                pure_tone['ix']=j
                            
                X = pure2cochlea( pure_tone, data )
                if func:
                    pure_tone['out']= func(X)
                else:
                    pure_tone['out']= X['Y']
                
                responses.append(pure_tone)

    if nthreads>len(f0_list):
        nthreads = len(f0_list)
        
    chunks_list = chunks( range(len(f0_list)),nthreads )
    
    thr = range(nthreads)
    for k in xrange(nthreads):
        jj = chunks_list[k]
        thr[k] = worker(jj)
        thr[k].start()
        
    for k in xrange(nthreads):   
        thr[k].join()
        
    return responses    
    
def pure2cochlea_explore_multiprocess(pure_tone,data,nprocs = 1,func=None):
   
    from multiprocessing import Process,Queue
    
    f0_list = pure_tone['f0']
    A_db_list = pure_tone['amplitude_db']
    duration_list =pure_tone['duration']
           
    def worker(jj,out_q):
        
        response = []
        
        for j in jj:
            pure_tone = {'amplitude_db':A_db_list[j],'f0':f0_list[j],'duration':duration_list[j]}
            pure_tone['ix']=j
            
            X = pure2cochlea( pure_tone, data )
            if func:
                pure_tone['out']= func(X)
            else:
                pure_tone['out']= X['Y']
            
            response.append(pure_tone)
        
        out_q.put(response)
                

    if nprocs>len(f0_list):
        nprocs = len(f0_list)
        
    chunks_list = chunks( range(len(f0_list)),nprocs )
    
    procs = []
    
    out_q = Queue()
   
    for k in xrange(nprocs):
        jj = chunks_list[k]
        p = Process(target=worker,args=(jj, out_q))
        
        procs.append(p)
        p.start()
   
    responses = []
    for k in range(nprocs):
        responses= responses + out_q.get()
        
    for p in procs:
        p.join()       

    return responses
    
def pure2periphery_explore_multiprocess(pure_tone,data,channels,nprocs = 1,func=None):
   
    from multiprocessing import Process,Queue
    
    f0_list = pure_tone['f0']
    A_db_list = pure_tone['amplitude_db']
    duration_list =pure_tone['duration']
           
    def worker(jj,out_q):
        
        response = []
        
        for i,j in enumerate(jj):
            pure_tone = {'amplitude_db':A_db_list[j],'f0':f0_list[j],'duration':duration_list[j]}
            pure_tone['ix']=j
            
            print i/float(len(jj))

            X = pure2periphery( pure_tone, data, channels )
            if func:
                pure_tone['out']= func(X)
            else:
                pure_tone['out']= X
            
            response.append(pure_tone)
        
        out_q.put(response)
                

    if nprocs>len(f0_list):
        nprocs = len(f0_list)
        
    chunks_list = chunks( range(len(f0_list)),nprocs )
    
    procs = []
    
    out_q = Queue()
   
    for k in xrange(nprocs):
        jj = chunks_list[k]
        p = Process(target=worker,args=(jj, out_q))
        
        procs.append(p)
        p.start()
   
    responses = []
    for k in range(nprocs):
        responses= responses + out_q.get()
        
    for p in procs:
        p.join()

    return responses

def synth2cochlea_explore_multiprocess(synth,synth_parameters,data,nprocs = 1,func=None):
   
    from multiprocessing import Process,Queue
    
    N = len(synth_parameters.values()[0])
               
    def worker(jj,out_q):
        
        response = []
        
        for j in jj:
            sp = {k:v[j] for k,v in synth_parameters.iteritems()}
            signal = synth(fs=data['fs'],**sp)
            X = runcochlea( signal, data )
            
            sp['ix']=j
            if func:
                sp['out']= func(X)
            else:
                sp['out']= X['Y']
            
            response.append(sp)
        
        out_q.put(response)
                

    if nprocs>N:
        nprocs = N
        
    chunks_list = chunks( range(N),nprocs )
    
    procs = []
    
    out_q = Queue()
   
    for k in xrange(nprocs):
        jj = chunks_list[k]
        p = Process(target=worker,args=(jj, out_q))
        
        procs.append(p)
        p.start()
   
    responses = []
    for k in range(nprocs):
        responses= responses + out_q.get()
        
    for p in procs:
        p.join()       

    return responses

def synth2periphery_explore_multiprocess(synth,synth_parameters,data,channels=0,nprocs = 1,func=None):
   
    from multiprocessing import Process,Queue
    
    N = len(synth_parameters.values()[0])
               
    def worker(jj,out_q):
        
        response = []
        
        for j in jj:
            sp = {k:v[j] for k,v in synth_parameters.iteritems()}
            signal = synth(fs=data['fs'],**sp)

            X = runcochlea( signal, data )

            Y = X['Y']
            syn,vihc = ihc_synapse(Y,data,channels=channels)

            out = {'Y':Y,'syn':syn}

            sp['ix']=j

            if func:
                sp['out']= func(out)
            else:
                sp['out']= out
            
            response.append(sp)
        
        out_q.put(response)
                

    if nprocs>N:
        nprocs = N
        
    chunks_list = chunks( range(N),nprocs )
    
    procs = []
    
    out_q = Queue()
   
    for k in xrange(nprocs):
        jj = chunks_list[k]
        p = Process(target=worker,args=(jj, out_q))
        
        procs.append(p)
        p.start()
   
    responses = []
    for k in range(nprocs):
        responses= responses + out_q.get()
        
    for p in procs:
        p.join()       

    return responses