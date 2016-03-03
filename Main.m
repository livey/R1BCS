%% This is script test the two algorithm R1BCS , 1BCS (oneBCS)
clear all;close all;clc;
% Macros 
M=100;
N=50;
K=5;
L=3;
max_iter=500;

[t,A,x,loc]=Gen_Data_Flip(M,N,K,L);
[x_est1,mun]=R1BCS(t,A,max_iter);
[x_est2]=oneBCS(t,A,max_iter,1e-6);

%% plot original signal and the estimated signal 
figure;
stem(x,'o','linewidth',2);
hold on;
stem(x_est1,'sr','linewidth',2);
stem(x_est2,'pc','linewidth',1.5);
h=legend('Original','R1BCS','1BCS');
set(h,'fontsize',10);

%% plot the sign flip location and estimated mun,which indicates the 
%   possible sign flip locations 
 temp=zeros(N,1);
 temp(loc)=1;
 figure;
 stem(temp,'o','linewidth',2);
 hold on;
 stem(mun,'rs','linewidth',2);
h=legend('sign-flip location','\mu_n');
set(h,'fontsize',12);
 
