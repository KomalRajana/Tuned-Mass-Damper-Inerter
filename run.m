% EXAMPLE
clc
clear
%%
m1=5; % Mass of primary
mTMD=0.25; %Mass of TMD
M1=[mTMD 0;0 m1] % Tonnes
K1=19800; % KN/m
omg1=(K1/m1)^0.5; %rad/s
C1=2*0.05*omg1*m1; % Arbitrary
B=2.5;
M=M1+[B 0;0 0]

% x=[0.1  1];

%%
% G1TMD=minmax_fun_SDOF(x,M, K1, C1,M1,B,omg1,mTMD);

%%

%% Optimization Algorithms
% common parameters
b=[]; %not Inertance
A=[];
Aeq=[];
beq=[];

funs=@(x)minmax_fun_SDOF(x,M, K1, C1,M1,B,omg1,mTMD);

x0=[0.01 0.01];
lb=[0.01 0.01]; % lower bound for stiffness and damping
ub=[1 1]; % upper bound for stiffness and damping

% Pattern Search % Plot Function== 'PlotFcn','psplotbestf',
options = optimoptions('patternsearch','UseParallel', false, 'UseCompletePoll', false, 'UseVectorized', false); % Minimize abs. values
[x,fval] =patternsearch(funs,x0,A,b,Aeq,beq,lb,ub,options);