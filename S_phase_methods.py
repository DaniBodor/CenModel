'''
During S phase, CENP-A is being redistributed at random over 2 nasCent Chrosomes
We will only follow one of these Chrosomes
Computationally, each CA(1) has a 50% chance of being replaced by an H3(0)
'''



#### what is faster? nprnd.rand() or nprnd.randint()??
import numpy as np
import numpy.random as nprnd
from timeit import timeit as timer 
from random import random, randrange



#creates a random chromosome of same size as original and does bitwise_and operation on the two
def Sphase_fullRandChro(Chro):
    randChro = nprnd.randint(2,size=ChroL)	# creates a random chromo of 0s and 1s. IS THERE A WAY TO ONLY WORK WITH THE 1s OF ORIGINAL CHROMO??
    newChro = np.bitwise_and(Chro,randChro)		# only returns a 1 if old Chro and randChro have a 1, so 50% of Chro-1s are replaced by 0
    
    return newChro

#creates a random 0-OR-1 for each CA and operates only on original CAs
def Sphase_onlyCAs(Chro):	#NOT WORKING!
    totCAs = Chro.count(1)
    randChro = nprnd.randint(2,size=totCAs)	# creates a random 0 or 1 for each existing CENP-A

    #ChroCAs = filter(lambda rm: rm != 0, Chro)
    newChro = [0 for CA in Chro if nprnd.randint(2)]

    return newChro

	
	
	
def Sphase_oneloop(Chro):
    n=0
    for CA in Chro:
        if CA == 1:
            if nprnd.rand() < 0.5:
                Chro[n] = 0
        n+=1
    return Chro

	
	
	
#creates a random chromosome of (0:1)-floats in a loop and multiplies (+round+int) each element with original chromo in a for loop
def Sphase_forloops(Chro):
    oldChro = Chro
    randChro = [randrange(0,1e10)/1e10 for _ in xrange(ChroL)]
    newChro = [int(round(Chro*randChro,0)) for Chro,randChro in zip(Chro,randChro)]

    return newChro
	