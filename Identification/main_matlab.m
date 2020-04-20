clear 
close all
clc

%% PLANT MODEL
load Identified3_real_inputs
H=B(:,3);
B=B(:,1:2);
D=D(:,1:2);
A = A(1:4,1:4);
B = B(1:4,:);
H = H(1:4,:);
C = C(:,1:4);

%%
Bw=H;
sys1=ss(A,Bw,C,0);

%%%%%%% Kalma Filter Design
sysk=ss(A,[B Bw],C,0);

[kest,L,P] = kalman(sysk,1e1,eye(2),0);

%%%%%%% Optimal Gain
syst=ss(A,B,C,0);
[K,S,e] = lqry(syst,1e2*eye(2),1e0*eye(2),0);

%%%%%%% Controller
rlqg = lqgreg(kest,K);
