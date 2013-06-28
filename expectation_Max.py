# Expectation maximization for Localization 
import math
import time
import random
import pickle


#Grid [11*7]
x_max=10
y_max=10

'''
Latent Variables:
1.Position of object (one of the grid points) Dimension J=11*7
2.Power of Transmitter (accounting for 10dBm of power variations with 0.5 dBm granularity) Dimension K=20


'''
J=x_max*y_max #dimension of position vector
K=20 #dimension of power levels
N=4 #number of sniffers 
powerLevelMin=20 #dBm
'''
Model Parameters:
X=probability distribution of being at any location in the grid
PL=probability distribution of having one of the power level 
Received power is modeled as gaussian  
mu[i][j][k] : mean of received power at ith sniffer from jth position and kth power level
sigma[i][j][k] : variance of received power at ith sniffer from jth position and kth power level
'''
X=[0 for i in range(J)]
PL=[0 for i in range(K)]
mu=[[[0 for k in range(K)] for j in range(J)] for i in range(N)]
sigma=[[[0 for k in range(K)] for j in range(J)] for i in range(N)]
PL_EXP=[0 for i in range(N)]
 
# initalization const
sniffer_loc=[(0,0),(0,y_max),(x_max,0),(x_max,y_max)]
d0=1
#

#EM details
MAX_ITER=5
# Dimension of measurements
M=10 
measurements=[[0 for i in range(N) ] for j in range(M)]
#

def init():
	#initializing X with uniform distribution
	for i in range(J):
		X[i]=(1.0/J)
	#initializing PL with uniform distribution
	for i in range(K):
		PL[i]=(1.0/K)
	#initializing sigma with constant 6dBm 
	for i in range(N):
		for j in range(J):
			for k in range(K):
				sigma[i][j][k]=6
	#initializing path loss exponent for each sniffer to 3 (experimental value for indoor environment)
	for i in range(N):
		PL_EXP[i]=3
	#initializing mu according to LDPL model
	for i in range(N):
		for j in range(J):
			for k in range(K):
				d=math.sqrt(((j%x_max)-sniffer_loc[i][0])**2+((j/x_max)-sniffer_loc[i][1])**2)
				if (d<d0):
					mu[i][j][k]=powerLevelMin
				else:
					mu[i][j][k]=powerLevelMin+0.5*k-10.0*PL_EXP[i]*math.log10(d/d0)

def Gaussian(m, s, x):
	# calculates the probability of x for 1-dim Gaussian with mean mu and var. sigma
	return math.exp(- ((m - x) ** 2) / (s ** 2) / 2.0) / math.sqrt(2.0 * math.pi *
	(s ** 2))

def probS_jk(j,k,l):
	acc=1
	for i in range(N):
		acc*=Gaussian(mu[i][j][k],sigma[i][j][k],measurements[l][i])
	return acc

def e_step():
	# Using current variable finding Expected values of latent variable
	theta=[[[0 for k in range(K)]for j in range(J)] for i in range(M)]

	
	for m in range(M):
		acc=0.0
		for j in range(J):
			for k in range(K):
				theta[m][j][k]=X[j]*PL[k]*probS_jk(j,k,m)
				acc+=theta[m][j][k]
		#normalizing
		for j in range(J):
			for k in range(K):
				theta[m][j][k]=theta[m][j][k]/acc
	return theta

def update_step(theta):
	#updating position estimate
	for j in range(J):
		acc=0.0
		for m in range(M):
			for k in range(K):
					acc+=theta[m][j][k]
		X[j]=acc/M
	
	#updating power level distribution
	
	for k in range(K):
		acc=0.0
		for m in range(M):
			for j in range(J):
				acc+=theta[m][j][k]
		PL[k]=acc/M

	# updating mu and sigma
	for i in range(N):
		for j in range(J):
			for k in range(K):
				acc_mu=0.0
				acc_sigma=0.0
				normalizing_const=0.0
				for m in range(M):
					acc_mu+=theta[m][j][k]*measurements[m][i]
					acc_sigma+=theta[m][j][k]*(measurements[m][i]-mu[i][j][k])*((measurements[m][i]-mu[i][j][k]))
					normalizing_const+=theta[m][j][k]
				#print normalizing_const
				mu[i][j][k]=acc_mu/normalizing_const
				#sigma[i][j][k]=max(math.sqrt(acc_sigma)/normalizing_const,0.5)

def printVariables(var_name):
	if('posD' in var_name):
		print "Position Distribution :"
		for i in range(y_max):
			for j in range(x_max):
				print X[i*x_max+j],
			print
	elif("powerD" in var_name):
		print "Power Distribution :"
		for i in range(K):
			print PL[i],
	elif("mus" in var_name):
		print "Mu's :"
		for i in range(N):
			for j in range(J):
				for k in range(K):
					print mu[i][j][k],
				print
			print
		print
	elif("sigmas" in var_name):
		print "sigma's :"
		for i in range(N):
			for j in range(J):
				for k in range(K):
					print sigma[i][j][k],
				print
			print
		print
	elif("data" in var_name):
		for i in range(M):
			for j in range(N):
				print measurements[i][j],
			print
		print
	else:
		pass


def generate_data(transmit_power,path_loss_exp,fading_var,loc):
	for i in range(M):
		for s in range(len(sniffer_loc)):
			d=math.sqrt((loc[0]-sniffer_loc[s][0])**2+(loc[1]-sniffer_loc[s][1])**2)
			if(d<d0):
				measurements[i][s]=transmit_power
			else:
				measurements[i][s]=transmit_power-10*path_loss_exp*math.log10(d/d0)+random.gauss(0,fading_var)
		

def main():
	NS=50
	error_exp=0
	error_round=0
	error_mprob=0
	for samples in range(NS):
		#initializing model parameters
		init()
		d_pos_x=random.randrange(x_max)
		d_pos_y=random.randrange(y_max)
		generate_data(22,4,6,(d_pos_x,d_pos_y))
		#printVariables(["data"])
		pos_x=0
		pos_y=0
		delta=1e-3
		it=0
		while(it<20):
			update_step(e_step())
			expected_pos_x=0.0
			expected_pos_y=0.0
			mep=0.0
			for i in range(J):
				expected_pos_y+=(i/x_max)*X[i]
				expected_pos_x+=(i%x_max)*X[i]
				if(X[i]>mep):
					mep=X[i]
					mi=i
			#print expected_pos_r,expected_pos_c,it
			if((pos_x-expected_pos_x)**2+(pos_y-expected_pos_y)**2<delta ):
				break
			it+=1
			pos_x=expected_pos_x
			pos_y=expected_pos_y
		print "Data Pos : ",d_pos_x,d_pos_y
		print "Expected Pos:",expected_pos_x,expected_pos_y
		print "Expected round off pos :",round(expected_pos_x),round(expected_pos_y)
		print "Most Prob Pos: ",mi%x_max,mi/x_max
		print 
		error_exp+=math.sqrt((d_pos_x-expected_pos_x)**2+(d_pos_y-expected_pos_y)**2)
		error_round+=math.sqrt((d_pos_x-round(expected_pos_x))**2+(d_pos_x-round(expected_pos_x))**2)
		error_mprob+=math.sqrt((d_pos_x-(mi%x_max))**2+(d_pos_x-mi/x_max)**2)
	print "MDE for Expected Pos: ",error_exp/NS
	print "MDE for Round off Pos: ",error_round/NS
	print "MDE for most prob Pos: ",error_mprob/NS
	#printVariables(["posD"])
	#print sum(X)



if __name__=="__main__":
		main()