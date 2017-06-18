clc
clear all
close all
warning off
dbstop if error

%% Initialize Data

 All permutations of phase factor B
M=4;
NN = 100;       % Number of iterations
N = 16;         % number of subcarriers
L = 4;
r=randint(NN,N,[0 M-1]);
% oversampling factor

%% Initialize the process
p=[1 -1 i -i]; % phase factor possible values
B=[];
U(cc)=16;
B=randsrc(U(cc),N,p)