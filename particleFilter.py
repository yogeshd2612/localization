#Particle filter implementation to track moving tx using RSSI measurements
from math import *
import random

world_size=(100,100)

#location of Sniffers at four corners
sniffer_loc=[(0,0),(0,100),(100,0),(100,100)]

#State of transmitter include its pos(x,y), power level, velocity in (x,y), 
#path_loss_exp,var due to evironment around each landmark/sniffer.


K=20
V_max=5
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
		t_m=mobile_tx()
		t_m.x=(self.x+self.v_x)%world_size[0]
		t_m.y=(self.y+self.v_y)%world_size[1]
		t_m.v_x=self.v_x
		t_m.v_y=self.v_y
		t_m.p=self.p
		t_m.pl_exp=self.pl_exp
		t_m.var=self.var
		return t_m
		
	def Gaussian(self,mu,sigma,x):
		# calculates the probability of x for 1-dim Gaussian with mean mu and var. sigma 
		return exp(- ((mu - x) ** 2) / (sigma ** 2) / 2.0) / sqrt(2.0 * pi * (sigma ** 2))

	def measurement_prob(self,measurements):
		prob=1.0
		X=self.sense()
		for i in range(len(sniffer_loc)):
			#prob*=self.Gaussian(X[i],self.var[i],measurements[i])
			prob*=self.Gaussian(X[i],5,measurements[i])
		return prob

	def __repr__(self):
		return '[x=%.6s y=%.6s p=%.6s]' % (str(self.x),str(self.y),str(self.p))


# Number of Particles
N=10000
# time units
T=20
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
	f=open("data.txt","w")
	f.write(str(T)+" "+str(N)+' 0 0 0 \n')
	for t in range(T):
		print "error : ", error(pos_data[t],p)
		'''print "Initial"
		for i in range(N):
			print p[i].x,p[i].y,p[i].v_x,p[i].v_y
		print''' 
		#writing data to file
		f.write(str(pos_data[t][0])+ " "+ str(pos_data[t][1])+" 0 0 0 \n")
		for i in range(N):
			f.write(str(p[i].x)+" "+str(p[i].y)+" "+str(p[i].p)+" "+str(p[i].v_x)+" "+str(p[i].v_y)+ '\n')
		
		# measurement step
		w=[0 for i in range(N)]
		for j in range(N):
			w[j]=p[j].measurement_prob(measurement_data[t])
		normalizer=sum(w)
		for i in range(N): w[i]/normalizer

		#print sorted(w)
		# resampling
		index = int(random.random()*N)
		beta =0.0
		mw=max(w)
		#print "max weight : ",mw
		p1=[]
		weight_avg=0.0
		for i in range(N):
			beta+=random.random()*2.0*mw
			while beta> w[index]:
				beta -= w[index]
				index=(index+1)%N
			weight_avg+=w[index]
			p1.append(p[index])
		'''print "after selection"
		for i in range(N):
			print p1[i].x,p1[i].y,p1[i].v_x,p1[i].v_y
		print''' 
		print "avg_weight of selection  : ",weight_avg/N
		# move step
		p2=[]
		for i in range(N):
			p2.append(p1[i].move())

		'''print "after move"	
		for i in range(N):
			print p2[i].x,p2[i].y,p2[i].v_x,p2[i].v_y
		print'''

		p=p2
	f.close()		

if __name__=="__main__":
	main()
