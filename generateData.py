import math
import random
import time
import pickle
#grid [11*7] granularity of 0.5

sniffer_loc=[(0,0),(3,0),(0,5),(3,5)]
transmit_power=25 #dBm at 0.5
d0=0.5
path_loss_exp=4
noise_sigma=6 #dbm
measurements=[]

def main():
	rows=11
	cols=7
	for i in range(rows):
		measurements.append([])
		for j in range(cols):
			ls=[]
			for s in range(len(sniffer_loc)):
				d=math.sqrt((j*0.5-sniffer_loc[s][0])**2+(i*0.5-sniffer_loc[s][1])**2)
				if(d<d0):
					ls.append(transmit_power)
				else:
					#Using LDPL model with some randomness for multipath effects
					ls.append(transmit_power-10*path_loss_exp*math.log10(d)+random.gauss(0,noise_sigma))  
			measurements[i].append(ls)

	#print measurements
	#Write to file
	f=open("data","w")
	pickle.dump(measurements,f)


if __name__=="__main__":
	main()