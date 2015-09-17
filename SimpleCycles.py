#Import libraries, functions, etc
#	and set parameters
import numpy as np
import numpy.random as nprnd
import pylab as plt
from S_phase_methods import Sphase_oneloop as Sphase2

ChroL = int(1000)					# length of chromosome (=Chro) in nucleosomes
initCA_Chro = 0.1					# fraction of initial CA on entire Chro
CenL = ChroL/100                		# UNUSED # length of centromere (=Cen) in nucleosomes
initCA_Cen = 0.3    				# UNUSED # fraction of initial CA on Cen
initCenPos = 0.5					# UNUSED # defines where Cen is located on Chro (from left to right)
erosion = 1e-2					# UNUSED # CA loss other than redistribution
divs = 25						# number of cycles through S and G1 phase
PlotC = 0                               # the plot counter ensures that separate plots are not overwritten


#%%


# define functions to run on Chrosome

# During S phase, CENP-A is being redistributed at random over 2 nasCent Chrosomes
# We will only follow one of these Chrosomes
# Computationally, each CA(1) has a 50% chance of being replaced by an H3(0)
def Sphase(Chro):
    newChro = Chro
    n=0
    for CA in newChro:
        if CA == 1:
            if nprnd.rand() < 0.5:
                newChro[n] = 0
        n+=1
    return newChro



# During G1phase, CENP-A is replenished onto empty (H3) nucleosome positions
# Rules for replenishment remain unclear!
#	need some density dependent rules to recruit CA to the proper site in the Chrosome
def G1phase(Chro):
    ''' THIS IS THE PLACE WHERE I NEED TO PLAY AROUND TO MAKE THE FUCKING MODEL WORK'''
    newChro = Chro
    fracCA = float(Chro.count(1))/len(Chro)
    #print fracCA
    n=0
    for CA in newChro:
        if CA == 0:
            if nprnd.rand() < 1.5*fracCA:
                newChro[n] = 1
        n+=1
    return newChro



# Create plot of CENP-A fraction
def FractionPlot(initialFraction, PlotC) :
    PlotC +=1
    plt.figure(PlotC)
    plt.title('Figure %i: positions occupied by CENP-A (out of %i)'%(PlotC,ChroL))

    xFP = np.arange(1+divs*2)
    yFP1 = ChroCounts
    yFP2 = [initialFraction]*(divs*2+1)
    plt.plot(xFP, yFP1, xFP, yFP2)
    
    plt.xlabel('Division')
    plt.ylabel('Fraction')
    my_xticks = range(divs+1)
    plt.xticks(xFP*2,my_xticks)
    plt.axis([0, 2*divs, 0, 1])

    plt.show()
    return PlotC
    

# Create graph of CENP-A positions on Centromere
def CenPlot(Chro,PlotC):
    PlotC +=1
    plt.figure(PlotC)    
    plt.title('Figure %i: chromosome usage by CENP-A'%(PlotC))
    
    xCP = np.arange(ChroL)
    yCP = Chro
    plt.bar(xCP, yCP, linewidth=0)
    
    plt.xlabel('Nucleosome position')
    plt.ylabel('Times occupied (in %i cycles)'%(divs))
#    my_xticks = range(divs+1)
#    plt.xticks(xFP*2,my_xticks)
    plt.axis([0, ChroL, 0, divs])
    
    plt.show()
    return PlotC
#%%


# Initialization of Chrosome. 0=H3 , 1=CENP-A
Chro = [0]*ChroL

for a in xrange(ChroL) :		# is there a way to replace this loop by a list comprehension?
    if nprnd.rand() < initCA_Chro :
        Chro[a] = 1
initialFraction = (float(Chro.count(1))/len(Chro))

# Initialization of Centromere
'''for b in xrange(CenL) :		# is there a way to replace this loop by a list comprehension?
    if nprnd.rand() < initCA_Chro :
        Chro[b+initCenPos/2] = 1
initialCenFraction = (float(Chro.count(1))/len(Chro))
'''


#%%

# looping through consecutive rounds of S- and G1phase

#print Chro#, Chro.count(1)

ChroCounts=[]
ChroCounts.append(float(Chro.count(1))/len(Chro))

CAatPos=np.array([0]*ChroL)

for a in xrange(divs) :
    Chro = Sphase2(Chro)
    #print Chro, Chro.count(1)
    ChroCounts.append(float(Chro.count(1))/len(Chro))
    if Chro.count(1) == 0: break
    
    Chro = G1phase(Chro)
    #print Chro, Chro.count(1)
    ChroCounts.append(float(Chro.count(1))/len(Chro))
    if Chro.count(1) == 0: break
    
    CAatPos += Chro

#%%
PlotC = FractionPlot(initialFraction,PlotC)
#%%
PlotC = CenPlot(CAatPos,PlotC)

#print Chro

#print newChro
#print Chro.count(1),newChro.count(1),round(float(newChro.count(1))/Chro.count(1),2)
#print conv_rate
#print round(sum(conv_rate)/len(conv_rate),3)