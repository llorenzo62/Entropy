#!/usr/bin/env python
# coding: utf-8

# In[24]:


import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'notebook')


def PE(serie,m,tau=1,okfr=False,oktie='order'):
    '''
    serie: time serie. list, np.array or compatible
    m: embedded dimension
    tau: embedded delay
    okfr: if True outpus pattern's empirical distribution as a dictionary
    oktie: tie treatment. 
        'equal' same rank for all elements
        'order' assign rank by order or aparition
        'drop' drops out patterns with ties
        
    '''
    freq={}
   
    serie=list(serie)
    
    #If delay > 1
    if tau>1:
        ser=[]
        for t in range(tau):
            ser+=serie[t::tau]
        serie=ser
    
    #Counting patterns
    if oktie=='order':
        for index in range(m,len(serie)+1):
            #encoding
            y=[(val,key) for key,val in enumerate(serie[index-m:index])]
            key=[-1]*len(y)
            for index,item in enumerate(sorted(y)):
                key[item[1]]=index
            key=tuple(key)
            #counting
            freq[key]=freq.get(key,0)+1
    else:
        for index in range(m,len(serie)+1):
            #encoding
            s={val:indx for indx,val in enumerate(sorted(serie[index-m:index]))}
            key=tuple(s[item] for item in serie[index-m:index])
            #counting
            freq[key]=freq.get(key,0)+1
            
    #Drop out malformed patterns
    freq={key:val for key,val in freq.items() if len(key)==m}
    
    if oktie=='drop': #drop out ties
        #l0=len(freq)
        freq={key:val for key,val in freq.items() if len(key)==len(set(key))}
        #dl=l0-len(freq)
        #print('{} ties removed, remain {} patterns'.format(dl,len(freq)))
    
    #Compute PE
    f=np.array(list(freq.values()))
    f=f/f.sum()
    if okfr:
        return (-f*np.log(f)).sum(), freq
    return (-f*np.log(f)).sum()   

