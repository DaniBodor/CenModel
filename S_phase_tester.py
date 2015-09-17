#Import libraries, functions, etc
#	and set parameters
import operator
from random			import random, sample, randrange

	#probably good idea to install numpy and figure that one out!

chromo_length = 25		# length of chromosome in nucleosomes
initCA = 0.3			# fraction of initial CA



# define functions to run on chromosome

# During S phase, CENP-A is being redistributed at random over 2 nascent chromosomes
# We will only follow one of these chromosomes
# Computationally, each CA(1) has a 50% chance of being replaced by an H3(0)

def Sphase(Chromo):		# phase in which CENP-A is redistributed
	oldChromo = Chromo
	randChromo = [randrange(0,1e10)/1e10 for _ in range(chromo_length)]		# Would prefer doing this without for loop!
	newChromo = [int(round(Chromo*randChromo,0)) for Chromo,randChromo in zip(Chromo,randChromo)]		# would prefer doing this list multiplication without using for loop; need numpy for this

	try:
		return float(newChromo.count(1))/oldChromo.count(1)
	except ZeroDivisionError:
		return -1

	
conv_rate=[]

for test in xrange(30000):

	# Initialization of chromosome. 0=H3 , 1=CENP-A
	Chromo = [0]*chromo_length
	for a in xrange(chromo_length) :
		if random() < initCA :
			Chromo[a] = 1
	#print(Chromo)
	#Sphase(Chromo)
	val = Sphase(Chromo)
	if val >= 0:
		conv_rate.append(val)
	



#print Chromo
#print newChromo
#print Chromo.count(1),newChromo.count(1),round(float(newChromo.count(1))/Chromo.count(1),2)
#print conv_rate
print round(sum(conv_rate)/len(conv_rate),3)