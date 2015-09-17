# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 16:53:13 2015

@author: Dani
"""


'''
# During G1phase, CENP-A is replenished onto empty (H3) nucleosome positions
# Rules for replenishment remain unclear!
#	need some density dependent rules to recruit CA to the proper site in the Chrosome
'''
import numpy as np
import pylab as plt

import numpy.random as nprnd
from timeit import timeit as timer 
from random import random, randrange
from scipy.ndimage.filters import gaussian_filter1d as gf1d
from scipy.stats import gaussian_kde as kernsmoo
import scipy.signal as sig
import sys
from scipy import stats






def RandomG1_woComprehension(Chro,rate,returnPDF=False,normalize=False):
    ''' Outputs a new Chromosome which has converted 0s to 1s in Chro ...
        ... at random with p=rate '''
    newChro = Chro
    fracCA = CountCA(Chro)
    n=0    
    for CA in newChro:
        if CA == 0:
            if nprnd.rand() < rate:
                newChro[n] = 1
        n+=1
    if returnPDF:
        return nprnd.random_sample(len(Chro))
    return newChro
    
    
def RandomG1(Chro,rate,returnPDF=False,normalize=False):    
    ''' Outputs a new Chromosome which has converted 0s to 1s in Chro ...
        ... at random with p=rate (untested with list comprehensio.... should work tho'''
    if returnPDF:
        return nprnd.random_sample(len(Chro))
    return [1 if nprnd.rand() < rate else nuc for nuc in Chro]



def Gauss(Chro,sigma,sigrange=0,returnPDF=False,normalize=False, mode='Sum'):
    if mode == 'Blur' :
        X = 'X'
    




def GaussBlurG1(Chro,sigma,sigrange=0,returnPDF=False,normalize=False):
    ''' Outputs a new Chromosome which has converted 0s to 1s in Chro ...
        ... with odds based on a gaussian blur of Chro'''
    blurred = gf1d( [float(i) for i in Chro] , sigma, mode='constant')  # list comprehension needed to convert ints to floats    
    if normalize:
        NormalizePDF(blurred,sum(Chro))
    if returnPDF:
        blurred = [max(0,i) for i in (blurred - Chro)]    # converts 1s in Chro to 0s in blurred #### because Chro is binary, it either doesn't affect (if 0) or pushes it to 0 (if 1 which is always bigger than blurred)
        if sigrange:
#            raise NameError('CHECK FIRST WHETHER CUTOFF METHOD IS CORRECT, IT IS CURRENTLY UNTESTED')
            cutoff = stats.norm.pdf(sigma*sigrange,scale=sigma)     ## I THINK THIS IS CORRECT, UNTESTED THOUGH
            blurred = [0 if i < cutoff else i for i in blurred]     # turn low values to 0 to decrease memory used in numpy array
        return blurred  
    return [1 if nprnd.rand() < nuc else 0 for nuc in Chro+blurred]
    


def GaussSumG1(Chro,sigma,sigrange=0,returnPDF=False,normalize=False):
    ''' Outputs a new Chromosome which has converted 0s to 1s in Chro ...
        ... with odds based on XXXXXXXXXXXXXXXXXXXXXX'''
    GaussArray = np.zeros( (sum(Chro),len(Chro)) )
    n=0
    for i in xrange(len(Chro)):
        if Chro[i] == 1:
            GaussArray[n] = getGauss(i,len(Chro),sigma,sigrange)
            n+=1
    gaussSum = np.sum(GaussArray, axis=0)
    if normalize:
        NormalizePDF(gaussSum,sum(Chro))
    if returnPDF:
        return gaussSum
    return [1 if nprnd.rand() < gaussSum[nuc] else nuc for nuc in Chro]


    
def GaussNegProdG1(Chro,sigma,sigrange=0,returnPDF=False,normalize=False):
    ''' Outputs a new Chromosome which has converted 0s to 1s in Chro ...
        ... with odds based on XXXXXXXXXXXXXXXXXXXXXX'''
    GaussArray = np.zeros( (sum(Chro),len(Chro)) )
    n=0
    for i in xrange(len(Chro)):
        if Chro[i] == 1:
            GaussArray[n] = getGauss(i,len(Chro),sigma,sigrange)
            n+=1
    gaussNegProd = NegProd(GaussArray)
    if normalize:
        NormalizePDF(gaussNegProd,sum(Chro))
    if returnPDF:
        return gaussNegProd
    return [1 if nprnd.rand() < gaussNegProd[nuc] else nuc for nuc in Chro]



#==============================================================================
#==============================================================================
# # This part contains functions that are called as part of G1 functions above
#==============================================================================
#==============================================================================

def getGauss(pos,ChroL,sigma,sigrange):
    ''' create an array with a gaussian surrounding pos and pos=0 '''
    singlePDF = np.zeros(ChroL)       # 0s left untreated in np arrays
    size = int(sigma*sigrange)
    if 0 < size < ChroL :
        if size % 2 == 0:  # i.e. is even
            size += 1     #gaussian becomes asymmetrical if even
    else:
        raise NameError('HAVENT IMPLEMENTED YET THAT SIGRANGE=0 or > ChroL LEADS TO UNBOUNDED GAUSS; this could also lead to asymmetry problems in even-sized Chro')
    gauss=sig.gaussian(size,sigma)
    gauss[size/2]=0
    gauss /= sum(gauss)

    #the following leads to end-effect problems near Chro borders!
    start = pos-(size/2)
    end = pos+(size/2)+1
    try:        # checks for end effect problems
        singlePDF [start:end] += gauss
    except ValueError:        # fixes end effect problems
        gauss,start,end = BorderFixer(singlePDF,gauss,start,end)
        singlePDF [start:end] += gauss
    return singlePDF

def BorderFixer(P,G,A,Z):
    ''' crops start and end of gauss so that it fits within in Chro if pos is too close to border
    #### NOTE: the sum of gauss is no longer 1, meaning there is some erosion near the boundaries,
    #### which is probably more realistic than renormalizing to 1 and increasing the effect of boundary CAs
    #### compared to non-boundary CAs'''
    if Z > len(P):
        G = G[:(len(P)-A)]
        Z = len(P)
    if A < 0:
        G = G[-A:]
        A = 0
    return G, A, Z



def NegProd(matrix):
    ''' outputs a 1d array of chance that each pos is drawn at least once if there 
    is an independent draw for each of the individual positions.'''
    NotMatrix = 1 - matrix          # If matrix is an array of chances, NotMatrix is the chance that it doesn't happen
    prodPDF = 1 - NotMatrix.prod(axis=0)    
    return prodPDF
    


def NormalizePDF(PDF,preCA):
    ''' Normalizes the imput PDF to sum to the number of pre-existing CAs'''
    PDF /= sum(PDF)
    PDF *= preCA
    return PDF
    
    
