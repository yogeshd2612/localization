%%initialization
clear all; close all; clc

fprintf('reading data ...\n')
data=load('data.txt');
T=data(1,1)
N=data(1,2)
axis([0,100,0,100],"manual")
for i=1:T
	px=data(2+(i-1)*(N+1),1);
	py=data(2+(i-1)*(N+1),2);
	X=data(3+(i-1)*(N+1):3+(i-1)*(N+1)+N-1,1);
	Y=data(3+(i-1)*(N+1):3+(i-1)*(N+1)+N-1,2);
	hold on;
	plot(X,Y,'r*','MarkerSize',10);
	plot(px,py,'bo','MarkerSize',10);
	pause(2);
	clf
end