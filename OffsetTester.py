# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 15:54:21 2015

@author: Dani
"""

import numpy as np
import numpy.random as nprnd
import pylab as plt
#from S_phase_methods import Sphase_oneloop as SphaseExt
from G1_phase_methods import *
import time, sys
#from matplotlib import rc
#rc('text', usetex=True)
from scipy.ndimage.filters import gaussian_filter1d as gf1d


#	Set parameters
ChroL = 200					# length of chromosome (=Chro) in nucleosomes; TRUE NUMBER ~6.5e5
initCA_Chro = 0.001   				# fraction of initial CA on entire Chro; TRUE NUMBER ~1e-3
CenL = int(5e3)                		# length of centromere (=Cen) in nucleosomes; TRUE NUMBER ~5e3
initCA_Cen = 0.04   				# fraction of initial CA on Cen; TRUE NUMBER ~0.04
initCenPos = 0.3    				# defines where Cen is located on Chro (from left to right)
erosion = 0     #or: initCA_Chro		# UNUSED # CA loss other than redistribution. Can also be introduced as 'temperature' to erode and reload at low frequency throughout
divs = 200	         				# number of cycles through S and G1 phase
critical = initCA_Chro                  # the level of CA at which the sim breaks because CEN is lost
sigma = 35                              # used for gaussian blur (or other gaussian function)
gausspow = 1                           # power to raise the gaussian to
gaussOffset = 50                         # used to offset peak of nascent CA recruitment, must be larger than sigma to displace max pos
start = time.time()



CApool = 3      # soluble CENP-A, which is realistically in the prox of 26% of 91000 (not necessarily in G1 though...)
efficiency = .5                                    # chance of each sol CA being incorporated; MUST BE NON-ZERO




def G1phase(Chro,G1method=GaussBlurG1,power=gausspow):
    ''' Gaussian blur (or other method) of Chro is made which functions as the blueprint for the chance of nascent
    CENP-A incorporation. 
    The Distribution is only made once at the beginning of G1 phase, meaning that nascent CENP-A does not 
    guide more nascent CENP-A incorporation. This makes the code a lot faster.
    A method was built in to ensure that there is no 'reloading' of CENP-A at any given position'''
    global sendChro    
    if gaussOffset:
        sendChro = np.zeros(ChroL)
        for nuc in xrange(ChroL):
            if Chro[nuc] == 1:
                if nuc >= gaussOffset:     sendChro[nuc-gaussOffset] += 0.5
                if nuc < ChroL-gaussOffset:    sendChro[nuc+gaussOffset] += 0.5
    else: sendChro = Chro
    
    Distrib = gf1d( [float(i) for i in sendChro] , sigma, mode='constant')
    Distrib = [max(0,i) for i in (Distrib - Chro)]
    
    
    if power != 1:                          # >1 makes the distribution favor clusters (at least in theory)
        Distrib = np.power(Distrib,[power]*len(Distrib))
    print sum(Distrib)
    Distrib /= sum(Distrib)/efficiency      # normailizes the distribution to sum to efficiency (i.e. random numbers above efficiency will not lead to CA incorporation)
    Cumulative = np.cumsum(Distrib)         # converts the distribution to a cumulative distribution list
    plt.plot(Distrib)
    print sum(Distrib)
    
    reloads = 0
    for n in xrange(CApool):                # for each molecule in the available pool of CA
        X = nprnd.rand()                    
        if X < efficiency:                  # test if it will be loaded or not
            newPos = np.searchsorted(Cumulative,X)
            while Chro[newPos]:             # ensures that previously loaded sites are not 'reloaded'; this also alleviates the necessity of converting 1s to 0s in GaussBlur
                X = nprnd.rand()/2
                newPos = np.searchsorted(Cumulative,X)
                reloads += 1                # counter to check how often a previously loaded site revisited
            Chro[newPos] = 1                # convert a 0 from the cumulative distribution to 1
    return Chro , reloads



Chro = np.zeros(ChroL)
Chro[100] = 1

#print 'Chro', Chro
newChro, R = G1phase(Chro)


#print 'sendChro', sendChro
