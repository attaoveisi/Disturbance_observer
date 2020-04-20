clear
clc

loadmatfile('A.mat');
loadmatfile('B.mat');
loadmatfile('H.mat');
loadmatfile('C.mat');
loadmatfile('D.mat');
A=A(1:4,1:4);
B=B(1:4,:);
C=C(:,1:4);
H=H(1:4,:);
//D=D(:,1);

m=rank(B);
n=size(A);
n=n(1,1);
qw=rank(H);
p=rank(C);


[Q,Sigmaw,R]=svd(H);
Sigma1=Sigmaw(1:qw,:);
Q1=Q(:,1:qw);
Q2=Q(:,qw+1:$);

Ab=Q'*A*Q;
Bb=Q'*B;
Ab11=Ab(1:qw,1:qw);
Ab12=Ab(1:qw,qw+1:$);
Ab21=Ab(qw+1:$,1:qw);
Ab22=Ab(qw+1:$,qw+1:$);
Bb1=Bb(1:qw,:);
Bb2=Bb(qw+1:$,:);
Cb=C*Q;
Cb1=Cb(:,1:qw);
Cb2=Cb(:,qw+1:$);
Cb1p=(Cb1'*Cb1)^-1*Cb1';
dum1=size(Cb1p);
dum2=dum1(1,2);
dum3=dum1(1,1);
dum4=dum2-dum3;
//if dum4 ~= 0
//    M=rand(dum4,dum2);
//    while rank([Cb1p;M]) ~= dum2;
        M=rand(dum4,dum2);
//    end
//    N=[Cb1p;M];
//else
//    N=[Cb1p];
//    M=eye(p,p);
//end
//loadmatfile('M.mat');
N=[Cb1p;M];

At2=Ab22-Ab21*Cb1p*Cb2;
Bt2=Bb2;
Gt2=Ab21*Cb1p;
Ct=M*(eye(p,p)-Cb1*Cb1p)*Cb2;
Ht=M*(eye(p,p)-Cb1*Cb1p);
dum5=size(Cb1p*Cb2);
dum6=dum5(1,2);
Qbx=Q*[Cb1p*Cb2;eye(dum6,dum6)];
dum7=size(Cb1p);
dum8=dum7(1,2);
dum9=size(Q);
dum10=dum9(1,2);
dum11=dum10-dum7(1,1);
Qby=Q*[-Cb1p;zeros(dum11,dum8)];
dum15=size(M);

//gama0=0.1;
//gama=10^8;

function [LME, LMI, OBJ]=DOC(XLIST)
[P1,P2,Ph2,Kh,Lh,gama]= XLIST(:)
LME=list(P1-P1',P2-P2',P2*B-B*Ph2)
LMI=list(-([A'*P2+P2*A+C'*C+B*Kh+Kh'*B',-B*Kh,P2*H,P2*B;-Kh'*B',Q2*At2'*P1*Q2'+Q2*P1*At2*Q2'-Q2*Ct'*Lh'*Q2'-Q2*Lh*Ct*Q2'+Q2*Q2',zeros(n,qw+m);H'*P2,zeros(qw,n),-gama*eye(qw,qw),zeros(qw,m);B'*P2,zeros(m,n+qw+m)]),P1,P2,gama)
OBJ=[]
endfunction

P1_0=1*eye(n-qw,n-qw);
P2_0=eye(n,n);
Ph2_0=eye(m,m);
Kh_0=rand(m,n);
Lh_0=rand(n-qw,dum15(1,1));
gama_0=1e1;

Init_guess=list(P1_0,P2_0,Ph2_0,Kh_0,Lh_0,gama_0);
//Init_guess=list(P1_0,P2_0,Ph2_0,Kh_0,Lh_0,eps1_0,eps2_0,gama_0);
Mbound=1;
abstol=0.3;
nu=1;
maxiters=200;
reltol=0.01;
//Ans_LMI=lmisolver(Init_guess,DOC,[Mbound,abstol,nu,maxiters,reltol]);
Ans_LMI=lmisolver(Init_guess,DOC);

P1=Ans_LMI(1);
P2=Ans_LMI(2);
Ph2=Ans_LMI(3);
Kh=Ans_LMI(4);
Lh=Ans_LMI(5);
gama=Ans_LMI(6);

K=Ph2^-1*Kh;
L=P1^-1*Lh;

Sai1=Qbx*(At2-L*Ct)-A*Qbx;
Sai2=Qbx*(Gt2+L*Ht)-A*Qby;
Sai3=Qby;
Sai4=Qbx*Bt2-B;
Sai5=Qbx*Q2';

savematfile('K.mat','K');
savematfile('L.mat','L');
savematfile('M.mat','M');
