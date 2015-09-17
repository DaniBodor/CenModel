# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 16:55:18 2015

@author: Dani
"""


def Initialize()
    # Initialization of Chrosome. 0=H3 , 1=CENP-A
    Chro = [0]*ChroL
    
    for a in xrange(ChroL) :		# is there a way to replace this loop by a list comprehension?
        if nprnd.rand() < initCA_Chro :
            Chro[a] = 1
    initialFraction = (float(Chro.count(1))/len(Chro))
    
    return Chro

    # Initialization of Centromere
    '''for b in xrange(CenL) :		# is there a way to replace this loop by a list comprehension?
        if nprnd.rand() < initCA_Chro :
            Chro[b+initCenPos/2] = 1
    initialCenFraction = (float(Chro.count(1))/len(Chro))