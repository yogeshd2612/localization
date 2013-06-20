# MAP estimation for Localization
import math
import time
import random
import pickle


#Grid [11*7]
rows=11
cols=7

'''
Latent Variables:
1.Position of object (one of the grid points) Dimension J=11*7
2.Power of Transmitter (accounting for 10dBm of power variations with 0.5 dBm granularity) Dimension K=20
Actual Observed Data:
S (received power at four sniffer located at corner of grid): Dimension N=4

'''
J=77 #dimension of position vector
K=20 #dimension of power levels
N=4 #number of sniffers 
S=[]
LOC=[]
powerLevelMin=20 #dBm
'''
Model Parameters:
X=probability distribution of being at any location in the grid
PL=probability distribution of having one of the power level 
Received power is modeled as gaussian  
mu[i][j][k] : mean of received power at ith sniffer from jth position and kth power level
sigma[i][j][k] : variance of received power at ith sniffer from jth position and kth power level
ploss_exp[i] : path loss exponent for the ith sniffer
'''
X=[]
PL=[]
mu=[]
sigma=[]
PL_EXP=[]
 
# initalization const
sniffer_loc=[(0,0),(3,0),(0,5),(3,5)]
d0=0.5
#

#EM details
MAX_ITER=2
#

def init():
	#initializing X with uniform distribution
	for i in range(J):
		X.append(1.0/J)
	#initializing PL with uniform distribution
	for i in range(K):
		PL.append(1.0/K)
	#initializing sigma with constant 6dBm 
	for i in range(N):
		sigma.append([])
		for j in range(J):
			sigma[i].append([])
			for k in range(K):
				sigma[i][j].append(6)
	#initializing path loss exponent for each sniffer to 3 (experimental value for indoor environment)
	for i in range(N):
		PL_EXP.append(3)
	#initializing mu according to LDPL model
	for i in range(N):
		mu.append([])
		for j in range(J):
			mu[i].append([])
			for k in range(K):
				d=math.sqrt(((j%cols)*0.5-sniffer_loc[i][0])**2+((j/cols)*0.5-sniffer_loc[i][1])**2)
				if d<d0:
					mu[i][j].append(powerLevelMin)
				else:
					mu[i][j].append(powerLevelMin+0.5*k-10*PL_EXP[i]*math.log10(d/d0))

def Gaussian(mu, sigma, x):
	# calculates the probability of x for 1-dim Gaussian with mean mu and var. sigma
	return math.exp(- ((mu - x) ** 2) / (sigma ** 2) / 2.0) / math.sqrt(2.0 * math.pi *
	(sigma ** 2))

def probS_jk(j,k):
	acc=1
	for i in range(N):
		acc*=Gaussian(mu[i][j][k],sigma[i][j][k],S[i])
	return acc

def e_step():
	# Using current variable finding Expected values of latent variable
	theta=[]
	for j in range(J):
		theta.append([])
		for k in range(K):
			theta[j].append(0)

	acc=0.0
	for j in range(J):
		for k in range(K):
			theta[j][k]=X[j]*PL[k]*probS_jk(j,k)
			acc+=theta[j][k]
	#normalizing
	for j in range(J):
		for k in range(K):
			theta[j][k]=theta[j][k]/acc
	return theta

def update_step(theta):
	#updating position estimate
	expected_lI=0
	max_ep=0
	for j in range(J):
		acc=0.0
		for k in range(K):
				acc+=theta[j][k]
		X[j]=acc
		if(X[j]>max_ep):
			max_ep=X[j]
			expected_lI=j

	#updating power level distribution
	expected_pI=0
	max_ep=0
	for k in range(K):
		acc=0.0
		for j in range(J):
			acc+=theta[j][k]
		PL[k]=acc
		if(PL[k]>max_ep):
			max_ep=PL[k]
			expected_pI=k

	expected_power=expected_pI*0.5+powerLevelMin
	expected_xpos=expected_lI%cols
	expected_ypos=expected_lI/cols
	#print ((expected_ypos,expected_xpos,expected_power))
	#updating path loss exponent 
	for i in range(N):
		d=math.sqrt(math.sqrt((expected_xpos*0.5-sniffer_loc[i][0])**2+(expected_ypos*0.5-sniffer_loc[i][1])**2))
		if(d<d0):
			continue
		PL_EXP[i]=(expected_power-S[i])/(10*math.log10(d/d0))

	#updating mu with new path loss exponent

	for i in range(N):
		for j in range(J):
			for k in range(K):
				d=math.sqrt(((j%cols)*0.5-sniffer_loc[i][0])**2+((j/cols)*0.5-sniffer_loc[i][1])**2)
				if d<d0:
					mu[i][j][k]=expected_power
				else:
					mu[i][j].append(expected_power-10*PL_EXP[i]*math.log10(d/d0))
	return ((expected_ypos,expected_xpos,expected_power))
def printVariables():
	print "Position Distribution :"
	for i in range(rows):
		for j in range(cols):
			print X[i*cols+j],
		print

	print "Power Distribution :"
	for i in range(K):
		print PL[i],
	print "\nPath Loss Exponent for each sniffer :"
	for i in range(N):
		print PL_EXP[i]
	'''for i in range(N):
		for j in range(J):
			for k in range(K):
				print mu[i][j][k],
			print
		print
	print'''


def main():
	
	
	#initializing model parameters
	init()
	#Reading measurement data
	f=open("data","r")
	measurements=pickle.load(f)
	y=random.randint(0,rows-1)
	x=random.randint(0,cols-1)
	globals()['S']=measurements[y][x]
	print "Actual Measurements : ",S
	print "Actual Position : ",y,",",x
	print "Estimates :"
	#printVariables()
	for it in range(MAX_ITER):
		print update_step(e_step())
		printVariables()
		


	#print measurements


if __name__=="__main__":
		main()