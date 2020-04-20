clear
clc

loadmatfile('A.mat');
loadmatfile('B.mat');
loadmatfile('H.mat');
loadmatfile('C.mat');
loadmatfile('observer.mat');

m=rank(B);
n=size(A);
n=n(1,1);
qw=rank(H);
p=rank(C);

function [LME,LMI,OBJ]=DRC(XLIST)
[delta2,delta3,deltah3,betah,Thx,Thw,gama6,gama7]= XLIST(:)

LME=list(delta2-delta2',delta3-delta3',delta3*B-B*deltah3)

LMI=list(-([-H'*betah'-betah*H,betah*(Gama-A-Phi*C),-2*betah*theta*C*A-Thw'*B',zeros(qw,qw),delta2;(betah*(Gama-A-Phi*C))',zeros(n,n),-Thx'*B',zeros(n,qw+qw);(-2*betah*theta*C*A-Thw'*B')',-B*Thx,C'*C+delta3*A+A'*delta3+B*Thx+Thx'*B',delta3*H+B*Thw,zeros(n,qw);zeros(qw,qw+n),(delta3*H+B*Thw)',-gama6*eye(qw,qw),zeros(qw,qw);delta2',zeros(qw,n+n+qw),-gama7*eye(qw,qw)]),delta2,delta3,gama6,gama7,-1-betah*H)

OBJ=[]

endfunction

delta2_0=eye(qw,qw);
delta3_0=eye(n,n);
deltah3_0=zeros(m,m);
betah_0=zeros(qw,n);
Thx_0=eye(m,n);
Thw_0=eye(m,qw);
gama6_0=1e1;
gama7_0=1e1;

Init_guess=list(delta2_0,delta3_0,deltah3_0,betah_0,Thx_0,Thw_0,gama6_0,gama7_0);

Mbound=100;
abstol=1e-3;
nu=1;
maxiters=500;
reltol=1e-3;

Ans_LMI=lmisolver(Init_guess,DRC,[Mbound,abstol,nu,maxiters,reltol]);
//Ans_LMI=lmisolver(Init_guess,DRC);

delta2=Ans_LMI(1);
delta3=Ans_LMI(2);
deltah3=Ans_LMI(3);
betah=Ans_LMI(4);
Thx=Ans_LMI(5);
Thw=Ans_LMI(6);
gama6=Ans_LMI(7);
gama7=Ans_LMI(8);

Tx=inv(deltah3)*Thx;
Tw=inv(deltah3)*Thw;
beta=inv(delta2)*betah;

savematfile('Tx.mat','Tx');
savematfile('Tw.mat','Tw');
savematfile('beta.mat','beta');

disp(spec(A+B*Tx))

