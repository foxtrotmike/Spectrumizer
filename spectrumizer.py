# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 22:00:20 2015
Computes the position dependent k-spectum feature representation in sparse format
Given a protein it can be used to obtain a windowed position dependent k-mer composition

Created by: Dr. Fayyaz Minhas, afsar <at> pieas [dot] edu dot pk

@author: fayyaz
"""
from random import choice
from itertools import product
import numpy as np
from scipy.sparse import lil_matrix
AA='ACDEFGHIKLMNPQRSTVWY'
class PDSpectrumizer:
    """
    Position Dependent Spectrumizer class
    """
    def __init__(self,k = 1, hws = 20, norm = True, fall = False):
        """
        k: [1] the value of k-mer
        hws: Half window size (The full window is of size 2*hws+1)
        The first window starts at the first character of the sequence.
        norm: Whether to normalize the feature vector of each window to have 
            max. norm of 1.0
        fall: Whether to allow "fall-off" (True) or not (False)
            i.e., compute the features of window at the end of sequence which
            have may not be full formed windows or not
        """
        self.k = k
        d=len(AA)**k
        self.hws = hws
        self.fws = 2*hws+1
        self.sdidx=dict(zip([''.join(a) for a in list(product(*([AA]*k)))],range(d)))  
        self.isdidx = None
        self.win  = range(2*hws+1)
        self.fall = fall
        if norm:
            self.v = 1/np.sqrt(2*hws+1)
        else:
            self.v = 1.0
    def spectrumize(self,s):
        """
        s: input protein sequence (must not contain any non-canonical characters)
        """
        D = self.sdidx
        hws = self.hws
        fws = self.fws
        k = self.k
        win = self.win
        d=len(AA)**k
        x = [D[s[i:i+k]] if s[i:i+k] in D else None for i in range(len(s)-k+1)]
        if self.fall:        
            odim = len(x)
        else:
            odim = int(len(x)-2*hws)
        f = lil_matrix((odim,d*(fws)))
        for i in range(odim):    
            loc = x[i:i+fws]
            iv = np.ravel_multi_index((win[:len(loc)],loc),(fws,d))
            f[i,iv] = self.v
        
        return f.tocsr()
    def __verify___(self,v):
        """
        Sanity checking -- used to convert back to sequence from the feature repn
        Can be used to construct motifs
        """        
        hws = self.hws
        k = self.k
        d=len(AA)**k
        if self.isdidx is None:
            self.isdidx = dict(zip(self.sdidx.values(),self.sdidx.keys()))
        isdidx = self.isdidx
        iloc,ik = np.unravel_index(v, (2*hws+1,d))
        return iloc, isdidx[ik]
        
if __name__=='__main__':
    s = ''.join([choice(AA) for _ in range(101)])   #random protein sequence 
    S = PDSpectrumizer(k = 1, hws = 10)
    f = S.spectrumize(s) #feature representation
    print s
    w = 0 #Let's try to re-construct the last complete window from the feature repn
    print ''.join([S.__verify___(v)[1] for v in f[w].nonzero()[1]]) # tester
    
   
   
    
    
