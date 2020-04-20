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
[delta3,deltah3,Thx,Thw,gama6]= XLIST(:)

LME=list(delta3-delta3',delta3*B-B*deltah3)

LMI=list(-([C'*C+delta3*A+A'*delta3+B*Thx+Thx'*B',delta3*H+B*Thw;(delta3*H+B*Thw)',-gama6*eye(qw,qw)]),delta3,gama6)
//LMI=list(-([zeros(n,n),-Thx'*B',zeros(n,qw);-B*Thx,C'*C+delta3*A+A'*delta3+B*Thx+Thx'*B',delta3*H+B*Thw;zeros(qw,n),(delta3*H+B*Thw)',-gama6*eye(qw,qw)]),delta3,gama6)

OBJ=[gama6]

endfunction

delta3_0=eye(n,n);
deltah3_0=zeros(m,m);
Thx_0=rand(m,n);
Thw_0=rand(m,qw);
gama6_0=0.01;

Init_guess=list(delta3_0,deltah3_0,Thx_0,Thw_0,gama6_0);

Mbound=1000;
abstol=1e-6;
nu=1;
maxiters=500;
reltol=1e-6;

//Ans_LMI=lmisolver(Init_guess,DRC,[Mbound,abstol,nu,maxiters,reltol]);
Ans_LMI=lmisolver(Init_guess,DRC);

delta3=Ans_LMI(1);
deltah3=Ans_LMI(2);
Thx=Ans_LMI(3);
Thw=Ans_LMI(4);
gama6=Ans_LMI(5);

Tx=inv(deltah3)*Thx;
Tw=inv(deltah3)*Thw;

savematfile('Tx.mat','Tx');
savematfile('Tw.mat','Tw');

disp([spec(A),spec(A+B*Tx)])

