#Particle filter implementation to track moving tx using RSSI measurements
from math import *
import random

world_size=(100,100)

#location of Sniffers at four corners
sniffer_loc=[(0,0),(0,100),(100,0),(100,100)]

#State of transmitter include its pos(x,y), power level, velocity in (x,y), 
#path_loss_exp,var due to evironment around each landmark/sniffer.


K=20
V_max=10
PL_EXP_max=5
P_range=10
P_min=20
VAR_max=10
class mobile_tx:
	def __init__(self):
		self.x=random.random()*world_size[0]
		self.y=random.random()*world_size[1]
		self.v_x=random.random()*V_max
		self.v_y=random.random()*V_max
		self.p=random.random()*P_range+P_min
		self.pl_exp=random.random()*PL_EXP_max
		self.var=[random.random()*VAR_max for i in range(len(sniffer_loc))]

	def set(self,x,y,v_x,v_y,p,pl_exp,var):
		if(x>world_size[1] or y>world_size[0] or v_x> V_max or v_y>V_max or p>P_min+P_range):
			raise ValueError ,'state out of bound'
		self.x=x
		self.y=y
		self.v_x=v_x
		self.v_y=v_y
		self.p=p
		self.pl_exp=pl_exp
		self.var=var

	def sense(self):
		Z=[]
		for i in range(len(sniffer_loc)):
			dist=sqrt((self.x-sniffer_loc[i][0])**2+(self.y-sniffer_loc[i][1])**2)
			Z.append(self.p-10*self.pl_exp*log10(dist))
		return Z

	def move(self):
		self.x=(self.x+self.v_x)%world_size[0]
		self.y=(self.y+self.v_y)%world_size[1]

	def Gaussian(self,mu,sigma,x):
		# calculates the probability of x for 1-dim Gaussian with mean mu and var. sigma 
		return exp(- ((mu - x) ** 2) / (sigma ** 2) / 2.0) / sqrt(2.0 * pi * (sigma ** 2))

	def measurement_prob(self,measurements):
		prob=1.0
		X=self.sense()
		for i in range(len(sniffer_loc)):
			prob*=self.Gaussian(X[i],self.var[i],measurements[i])
		return prob

	def __repr__(self):
		return '[x=%.6s y=%.6s p=%.6s]' % (str(self.x),str(self.y),str(self.p))


# Number of Particles
N=1000
# time units
T=40
pos_data=[0 for i in range(T)]
measurement_data=[0 for i in range(T)]
def generate_data(p,pl_exp,v_x,v_y,x,y,var):
	for t in range(T):
		Z=[]
		for i in range(len(sniffer_loc)):
			dist=sqrt((x-sniffer_loc[i][0])**2+(y-sniffer_loc[i][1])**2)
			Z.append(p-10*pl_exp*log10(dist)+random.gauss(0.0,var[i]))
		measurement_data[t]=Z
		pos_data[t]=(x,y)
		x=(x+v_x)%world_size[0]
		y=(y+v_y)%world_size[1]

def error(r,p):
	s=0.0
	# calculation mean error
	for i in range(len(p)):
		dx =(p[i].x-r[0] + (world_size[0]/2.0))%world_size[0] -(world_size[0]/2.0)
		dy =(p[i].y-r[1]+ (world_size[1]/2.0))%world_size[1] - (world_size[1]/2.0)
		err=sqrt(dx*dx+dy*dy)
		s+=err
	return s/len(p) 

def main():
	generate_data(25,4,2,3,20,20,[4,5,3,4])
	#print measurement_data
	#print pos_data
	p=[mobile_tx() for i in range(N)]

	for t in range(T):
		print "error : ", error(pos_data[t],p)
		for i in range(5):
			print p[i]
		# measurement step
		w=[0 for i in range(N)]
		for j in range(N):
			w[j]=p[j].measurement_prob(measurement_data[t])
		#print sorted(w)
		# resampling
		index = int(random.random()*N)
		beta =0.0
		mw=max(w)
		p1=[]
		for i in range(N):
			beta+=random.random()*2.0*mw
			while beta> w[index]:
				beta -= w[index]
				index=(index+1)%N
			p1.append(p[index])
		# move step
		for i in range(N):
			p1[i].move()
		
		p=p1
		

if __name__=="__main__":
	main()
