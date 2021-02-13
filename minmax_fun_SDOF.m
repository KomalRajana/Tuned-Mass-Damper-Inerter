function G1TMD=minmax_fun_SDOF(x,M, K1, C1,M1,b, omg1,mTMD)
% Fminmax function using state space system representation
%{
M = Mass Matrix of primary structure with inerter and TMD
K1= Stiffness primary structure
C1= damping matrix of primary structure
M1= Mass Matrix of primary structure with TMD only
b=  Inertance
omg1= first natural frequency of primary structure
mTMD= mass of TMD
tmd_floor= location of TMD floor
spec_type= type of spectra
dof_min= DOF to be minimized
%}

csiTMD=x(1)
Freqratio=x(2)

kTMDI=Freqratio^2*omg1^2*(mTMD+b)
cTMDI=2*csiTMD*(mTMD+b)*Freqratio*omg1
% n1=size(M,1);
% Assembly of Stiffness Matrix with TMDI
K=(zeros(size(M,1)))

K(1,1)=kTMDI;
K(1,2)=-kTMDI;
K(2,1)=-kTMDI;
K(2,2)=kTMDI+K1;

% Assembly of Damping Matrix with TMDI
C=(zeros(size(M,1)));
                                         
C(1,1)=cTMDI;
C(1,2)=-cTMDI;
C(2,1)=-cTMDI;
C(2,2)=cTMDI+C1; % Damping matrix
%save('MM.mat','M','C','K')
% C
% K
n=size(M,1)-1
%% Transfer Function
% g1(s)=c0.(sI-A)^(-1) * B :: I dimention 2n+2
C0=zeros(2*n+2,1)'; % Is (2n+2) length Vector
C0(2)=1;

I=eye(2*n+2) % is the Identitiy matrix of size 2n+2
% Obtain A Matrix
% A = [ 0(n+1)      I(n+1)
%       -M^(-1)*K   -M^(-1)*K]
% B=[0(n+1)
%   I(n+1)]

ae1=zeros(n+1)
ae2=eye(n+1)
ae3=-(M\K); %A\B= inv(B)*A
ae4=-(M\C);
ae5=zeros(n+1,1);
ae6=ones(n+1,1);
A=[ae1 ae2;ae3 ae4];
B=[ae5;(-M1/M)*ae6];

%% White Noise
% whiteNoiseFactor=1;
% g2=@(w) abs(C0*(((1i*w).*I-A)\B))^2*whiteNoiseFactor;
% G1TMD=integral(g2,-inf,inf,'ArrayValued',true);

% 
% Lyapunov Solution
BBT=(2*pi).*(B*B');
lyap_sol=lyap(A,BBT); % Lyapunov Solution
% 
G1TMD=lyap_sol(2,2)

end