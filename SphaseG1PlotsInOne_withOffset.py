#Import libraries, functions, etc

'''
STUFF TO ADD IN THE FUTURE

- create a model where the gaussblur is offset by X positions in both directions


- think about efficiency differently!
    - Should be a function of nearby CAs, rather than a constant for the entire chromosome
    - I.e. Distrib should not be normalized against any fixed number


- create a system where CA has differential stability in- and outside of CEN


- gaussian blur does not need to create 0s out of 1s in Chro (this is fixed by reload checker in G1 phase)
    - test the number of reload events if the above is true


- parameter test situation:
    - criteria:
            - Cen maintenance (in fraction CA)
            - No 2y Cen appearance (how to test this?)
    - for working parameter set: 
        - create dictionary with key being tuple of parameter sets and value is the final (and initial??) Chro
    - for working AND not-working parameter sets
        - export them to a txt file with an OK/FAIL key (can later be used to create 2D/3D representation of working sets)


ultimately:
- create situation with 46 chromosomes of correct sizes, all at the same time and check stability of all 46
- check situation of multiple initialized CENs per chromosome
- check situation of sudden increase/decrease of CApool
- check situation of CENpos being telocentric (i.e. CenPos close to 0)
- what happens at higher/lower CApools, CenSizes, ChroLens, etc.
- check for 99% CI of CA halflife in MBoC paper and include this value as erosion

'''




import numpy as np
import numpy.random as nprnd
import pylab as plt
from S_phase_methods import Sphase_oneloop as SphaseExt
from G1_phase_methods import *
import time, sys
from scipy.ndimage.filters import gaussian_filter1d as gf1d
#from matplotlib import rc
#rc('text', usetex=True)


#	Set parameters
ChroL = int(1e5)					# length of chromosome (=Chro) in nucleosomes; TRUE NUMBER ~6.5e5
initCA_Chro = 0.001   				# fraction of initial CA on entire Chro; TRUE NUMBER ~1e-3
CenL = int(5e3)                		# length of centromere (=Cen) in nucleosomes; TRUE NUMBER ~5e3
initCA_Cen = 0.04   				# fraction of initial CA on Cen; TRUE NUMBER ~0.04
initCenPos = 0.3    				# defines where Cen is located on Chro (from left to right)
erosion = 0     #or: initCA_Chro		# UNUSED # CA loss other than redistribution. Can also be introduced as 'temperature' to erode and reload at low frequency throughout
divs = 120	         				# number of cycles through S and G1 phase
critical = initCA_Chro                  # the level of CA at which the sim breaks because CEN is lost
sigma = 3                              # used for gaussian blur (or other gaussian function); must be smaller than gaussOffset to displace max pos
gaussOffset = 10                         # used to offset peak of nascent CA recruitment, must be larger than sigma to displace max pos
gausspow = 3                            # power to raise the gaussian to (1 means no effect)
start = time.time()
resetChro = True                     # if False, Chro will be initialized with the final outcome of the precious run


CApool = int(round(initCA_Chro*ChroL + initCA_Cen*CenL,0))      # soluble CENP-A, which is realistically in the prox of 26% of 91000 (not necessarily in G1 though...)
efficiency = .5                                    # chance of each sol CA being incorporated; MUST BE NON-ZERO
# using pre-pool as sol-pool at 0.5 efficiency maintains initial values


# Find out from which to which nucleosome the Centromere reaches
CenStart = int(initCenPos*ChroL) - CenL/2        #
CenEnd = CenStart + CenL


def CountCA (Chro):
    '''counts the number of 1s (=CA) in Chro (must be array like)'''
    if type(Chro) is np.ndarray:
        fraction = float( np.count_nonzero(Chro) )/len(Chro)
    elif type(Chro) is list:
        fraction = float( Chro.count(1) )/len(Chro)
    else:
        sys.exit('Error in CountCA function, Chro is not an array or list')
    return fraction


def format_e(n, maxdec=-1, minconvert=1e3):
    ''' convert float/int to nicely formatted scientific notation string'''
    if n < minconvert: return n
    else:
        a = '%e' % n
        pre = a.split('e')[0].rstrip('0').rstrip('.')
        if 0 <= maxdec < len(a)-2 :
            pre = pre[0:maxdec+2].rstrip('.')
        post = 'e' + str(int(a.split('e+')[1]))
        return pre+post



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
            if nprnd.rand() < 0.5-erosion:
                newChro[n] = 0
        n+=1
    return newChro



# During G1phase, CENP-A is replenished onto empty (H3) nucleosome positions
# Rules for replenishment remain unclear!
#	need some density dependent rules to recruit CA to the proper site in the Chrosome
def G1phase(Chro):
    ''' Gaussian blur (or other method) of Chro is made which functions as the blueprint for the chance of nascent
    CENP-A incorporation. 
    The Distribution is only made once at the beginning of G1 phase, meaning that nascent CENP-A does not 
    guide more nascent CENP-A incorporation. This makes the code a lot faster.
    A method was built in to ensure that there is no 'reloading' of CENP-A at any given position'''

    # Create Offset array which will be used to calculate the frequency distribution of nascent CA loading    
    if gaussOffset:
        OffsetChro = np.zeros(ChroL)
        for nuc in xrange(ChroL):
            if Chro[nuc] == 1:
                if nuc >= gaussOffset:     OffsetChro[nuc-gaussOffset] += 0.5
                if nuc < ChroL-gaussOffset:    OffsetChro[nuc+gaussOffset] += 0.5
    else: OffsetChro = np.array(Chro)

    # Create frequency distribution based on gaussian blurred
    Distrib = gf1d( [float(i) for i in OffsetChro] , sigma, mode='constant')     #this line is no longer outsourced to G1_phase_Methods file!
    Distrib = [max(0,i) for i in (Distrib - Chro)]                              # convert 1s in Chro to 0s in Distrib so that there is not reloading of CA at old positions

    if gausspow != 1:                          # >1 makes the distribution favor clusters (at least in theory)
        Distrib = np.power(Distrib,[gausspow]*len(Distrib))
    Distrib /= sum(Distrib)/efficiency      # normailizes the distribution to sum to efficiency (i.e. random numbers above efficiency will not lead to CA incorporation)
    Cumulative = np.cumsum(Distrib)         # converts the distribution to a cumulative distribution list
    
    reloads=0                       # counter to check how often a previously loaded site revisited
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

    
def G1phaseRecalc(Chro,G1method=GaussBlurG1):
    ''' Alternative version of the one above where the Gaussian is recalculated each time a CENP-A is loaded
    i.e. nascent CENP-A promotes assembly more nascent CENP-A.
    This may or may not be more accurate, but is definitely WAAAYYYYY SLOWER!'''
    for n in xrange(CApool):
        X = nprnd.rand()
        if X < efficiency:                          # test if it will be loaded or not
            Distrib = G1method(Chro,sigma,sigrange=0,returnPDF=True)    # retrieves a gaussian conversion of Chro, based on specific G1 method from other file
            Distrib /= sum(Distrib)/efficiency      # normailizes the distribution to sum to efficiency (i.e. random numbers above efficiency will not lead to CA incorporation)
            Cumulative = np.cumsum(Distrib)         # converts the distribution to a cumulative distribution list
            newPos = np.searchsorted(Cumulative,X)  # check where random number falls on this cumulative list
            Chro[newPos] = 1                        # and convert that site to a 1
    return Chro



#%%
# define output formats (graphs)


def FractionPlot(ChroCA,CenCA):
    ''' This function creates a graph of CA occupancy over time throughout Chro and at Cen '''
    
    Title = 'Fraction of positions occupied by CENP-A'    
    plt.figure(Title)
    plt.title(Title)

    #plot Chro fraction (green): current (solid), initial (dotted), and average (dashed)
    plt.plot(ChroCA, 'g', label='Chro (%s sites)'%format_e(ChroL,1))
    plt.axhline(ChroCA[0], ls=':', c='g')
    plt.axhline(np.mean(ChroCA), ls='--', c='g')

    #plot Cen fraction (magenta): current (solid), initial (dotted), and average (dashed)
    plt.plot(CenCA, 'm', label='Cen (%s sites)'%format_e(CenL,1))
    plt.axhline(CenCA[0], ls=':', c='m')
    plt.axhline(np.mean(CenCA), ls='--', c='m')
    
#    plt.axhline(critical, ls='--', c='k', label="critical")
    
    plt.legend(frameon = False, ncol = 2, loc=9)
    plt.xlabel('Division')
    plt.ylabel('Fraction')
    plt.show()
    

def OccupancyPlot(fracs):
    ''' This function creates a graph of how often each position has been used throughout all divisions'''
    
    Title = 'Chromosome occupancy per site'
    plt.figure(Title)
    plt.title(Title)

    x = xrange(ChroL)
    y = fracs
    plt.axvspan(CenStart, CenEnd, color='b', alpha=0.5, lw=0)
    plt.bar(x, y, lw=0, width=1, color='black')
    
    plt.xlabel('Nucleosome position')
    plt.ylabel('Times occupied (in %i cycles)'%(divs))
    plt.axis([0, ChroL, 0, 1])
    plt.show()


def ChroPlot(initial,final):
    ''' This function creates a graph of initial and final Chromosomes'''
    
    Title = 'Chromosome usage'
    plt.figure(Title)
    plt.title(Title)
    
    plt.plot(initial, 'r', lw=.05, label='initial')
    plt.plot(final, 'b', lw=.05, label='final')
    plt.axvspan(CenStart, CenEnd, color='none', lw=1.3, ec='k')
    
    plt.ylabel('CA')
    plt.xlabel('site')
    #plt.legend()       #labels invisible due to extremely thin lines   
    plt.show()


def BinnedChroPlot(initial, final, binsize=100, CenShading=[CenStart,CenEnd], name='Chromosome'):
    ''' This function creates a graph of initial and final Chromosomes binned per (binsize) '''
    
    Title = '%s usage, binned per %i'%(name,binsize)
    plt.figure(Title)
    plt.title(Title)
    
    global initialHeights   
    
    initialSize = len(initial)/binsize
    finalSize = len(final)/binsize
    
    initialHeights = [sum(initial[x:x+binsize]) for x in xrange(0,len(initial),binsize)] # idem, but binned in bargraph[x:x+binsize]) for x in xrange(0,ChroL,binsize)] #create binned initial Chro
    finalHeights = [sum(final[x:x+binsize]) for x in xrange(0,len(final),binsize)]        #create binned Chro

    if CenShading:
        if CenShading == [CenStart,CenEnd]: CenShading = [CenStart/binsize,CenEnd/binsize]
        plt.axvspan(CenShading[0], CenShading[1], color='0.85')
    plt.bar(xrange(initialSize), initialHeights, lw=0, width=1, color='r', label='initial', alpha=0.4)
    plt.bar(xrange(finalSize), finalHeights, lw=0, width=1, color='b', label='final', alpha=0.4)
    
    plt.ylabel('sites occupied (of %i)'%(binsize))
    plt.xlabel('bin number (%i sites per bin)'%(binsize))
    plt.legend()
    plt.show()
    
    
#%%
    
    
# Initialization of Chrosome. 0=H3 , 1=CENP-A
    
if resetChro:
    Chro = np.array([0]*ChroL)
    Cen = Chro[CenStart:CenEnd]
    
    for a in xrange(ChroL) :		# initialize low density CA throughout chromosome
        if nprnd.rand() < initCA_Chro :
            Chro[a] = 1
    
    # Initialization of Centromere
    for b in xrange(CenL) :		# initialize higher levels if CA at CEN 
        if nprnd.rand() < initCA_Cen :
            Chro[ CenStart + b ] = 1
    
    initialChro = list(Chro)
    initialCen = list(Cen)





#%%
# looping through consecutive rounds of S- and G1phase

ChroCounts= [CountCA(Chro)]  # initialize array of CA on Chro at the end of each cycle (i.e. post-G1)
CenCounts = [CountCA(Cen)]  # idem for Cen
reloadCounts = []           # initialize array for number of times a 'reload' event happens in G1

CAatPos=np.array([0]*ChroL)


postG1 = 0
postS = 0

precycling = time.time()
for currentDiv in xrange(divs) :
    '''S PHASE --> LOSS OF 50% OF CA'''
    preS = time.time()
    Chro = SphaseExt(Chro)
    postS += time.time()-preS
#    ChroCounts.append(CountCA(Chro))
#    CenCounts.append(CountCA(Cen))
    if CountCA(Cen) <= critical:
        ChroCounts.append(CountCA(Chro))
        CenCounts.append(CountCA(Cen))
        break         # end sim if CA levels at CEN are below critical values
    
    '''G1 PHASE --> RELOADING OF CA'''
    preG1 = time.time()
    Chro,reloads = G1phase(Chro)
    postG1 += time.time()-preG1
    
    ChroCounts.append(CountCA(Chro))
    CenCounts.append(CountCA(Cen))
    reloadCounts.append(reloads)
    
    CAatPos += Chro
cyclingTime = time.time()-precycling


#%%
# ask for output



#ChroPlot(initialChro,Chro)      # ask for plot of Chromosome (initial and final)

BinnedChroPlot(initialChro,Chro,binsize=1000,name='Chromosome') # idem, but binned in bargraph
#BinnedChroPlot(initialChro,Chro,binsize=100,name='Chromosome') # idem, but binned in bargraph
#BinnedChroPlot(initialCen,Cen,binsize=100,CenShading=False,name='Centromere') # idem, but binned in bargraph
BinnedChroPlot(initialChro[CenStart-1000:CenEnd+1000],Chro[CenStart-1000:CenEnd+1000],binsize=100,name='Extended Centromere',CenShading=[10,60]) # idem, but binned in bargraph

#OccupancyPlot(XXXX) #don't remember what to feed into this but I think this: CAatPos/float(divs)

FractionPlot(ChroCounts,CenCounts)  # ask graph of CENP-A over time



print '*********************************************'

print 'Sigma:', sigma
print 'Offset:', gaussOffset
if gausspow != 1: print 'Gaussian raised to power:', gausspow
if currentDiv+1 < divs:
    divs = currentDiv
    print '****************** SUBCRITICAL LEVELS (%g) REACHED AFTER %i DIVISIONS'%(critical,currentDiv)    
print '*********************************************'
print 'total number of reloads in', divs, 'divisions:', sum(reloadCounts)
print 'average G1 phase runtime:', round(postG1/(divs),5), 'seconds'
print 'average S phase runtime:', round(postS/divs,5), 'seconds'
print 'time elapsed while cycling:', round(cyclingTime,1), 'seconds'
print 'total time elapsed:', round(time.time()-start,1), 'seconds'

