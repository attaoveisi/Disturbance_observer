clear
clc

loadmatfile('A.mat');
loadmatfile('B.mat');
loadmatfile('H.mat');
loadmatfile('C.mat');
H=[H,H];
C=eye(6,6);

m=rank(B);
n=size(A);
n=n(1,1);
[pokh,qw]=size(H);
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

function [LME,LMI,OBJ]=UIO(XLIST)
[N,J,L,P1,P2,P3,P4,gama1]= XLIST(:)
dum_var1=Cb1p*C*A;
dum_var2=Cb1p*Cb2*L*Ht*C;
LME=list(Cb1p*C*B-Cb1p*Cb2*J-Bb1,dum_var1-dum_var2,P1-P1',P2-P2',P3-P3',P4-P4')
LMI=list(-([N'-Ct'*L'+N-L*Ct+P1,Ct'*L'*Cb2'*Cb1p'-Cb2'*Cb1p',At2+Gt2*C*Q2-L*Ht*C*Q2+L*Ct-N,Gt2*C*Q1-L*Ht*C*Q1,Bb2-J,zeros(n-qw,qw);Cb1p*Cb2*L*Ct-Cb1p*Cb2,P4,Ab12-Cb1p*C*A*Q2+Cb1p*Cb2*N+Cb1p*Cb2*L*Ht*C*Q2-Cb1p*Cb2*L*Ct,Ab11,zeros(qw,m),Sigma1*R'-Cb1p*C*H;At2'+Q2'*C'*Gt2'-Q2'*C'*Ht'*L'+Ct'*L'-N',(Ab12-Cb1p*C*A*Q2+Cb1p*Cb2*N+Cb1p*Cb2*L*Ht*C*Q2-Cb1p*Cb2*L*Ct)',P2,zeros(n-qw,qw+m+qw);Q1'*C'*Gt2'-Q1'*C'*Ht'*L',zeros(qw,qw+n-qw),P3,zeros(qw,m+qw);Bb2'-J',zeros(m,qw+n-qw+qw+m+qw);zeros(qw,n-qw),R*Sigma1'-H'*C'*Cb1p',zeros(qw,n-qw+qw+m),-gama1]),P1,P2,P3,P4,gama1,-N)
OBJ=[]
endfunction

N_0=-eye(n-qw,n-qw);
J_0=zeros(n-qw,m);
L_0=zeros(n-qw,p-qw);
P1_0=eye(n-qw,n-qw);
P2_0=eye(n-qw,n-qw);
P3_0=eye(qw,qw);
P4_0=eye(qw,qw);
gama1_0=1e1;

Init_guess=list(N_0,J_0,L_0,P1_0,P2_0,P3_0,P4_0,gama1_0);

Mbound=1;
abstol=1e-5;
nu=5;
maxiters=500;
reltol=1e-5;

Ans_LMI=lmisolver(Init_guess,UIO,[Mbound,abstol,nu,maxiters,reltol]);
//Ans_LMI=lmisolver(Init_guess,UIO);

N=Ans_LMI(1);
J=Ans_LMI(2);
L=Ans_LMI(3);
P1=Ans_LMI(4);
P2=Ans_LMI(5);
P3=Ans_LMI(6);
P4=Ans_LMI(7);
gama1=Ans_LMI(8);

savematfile('N.mat','N');
savematfile('J.mat','J');
savematfile('L.mat','L');
savematfile('M.mat','M');

Omega1=R*(Sigma1^(-1))*(Cb1p*Cb2*L*Ct-Ab12-Cb1p*Cb2*N);
Omega2=-R*(Sigma1^(-1))*Ab11;
Omega3=R*(Sigma1^(-1))*Cb1p*C*H;

function [LME,LMI,OBJ]=DRC1(XLIST)
[R1,R3,R4,R6,Rh1,Khx,Khw,gama4]= XLIST(:)

P211=-Ab12'*R1*Cb1p*Cb2+Cb2'*Cb1p'*Ab11'*R1*Cb1p*Cb2-Q2'*Khx'*Bb1'*Cb1p*Cb2+Cb2'*Cb1p'*Q1'*Khx'*Bb1'*Cb1p*Cb2+Cb2'*Cb1p'*Omega2'*Khw'*Bb1'*Cb1p*Cb2-Omega1'*Khw'*Bb1'*Cb1p*Cb2-Cb2'*Cb1p'*R1*Ab12+Cb2'*Cb1p'*R1*Ab11*Cb1p*Cb2-Cb2'*Cb1p'*Bb1*Khx*Q2+Cb2'*Cb1p'*Bb1*Khx*Q1*Cb1p*Cb2+Cb2'*Cb1p'*Bb1*Khw*Omega2*Cb1p*Cb2-Cb2'*Cb1p'*Bb1*Khw*Omega1+R3+Cb2'*Cb1p'*R4*Cb1p*Cb2;
P212=-Cb2'*Cb1p'*Bb1*Khx*Q1*Cb1p*Cb2+Cb2'*Cb1p'*Bb1*Khx*Q2+Cb2'*Cb1p'*Bb1*Khw*Omega1-Cb2'*Cb1p'*Bb1*Khw*Omega2*Cb1p*Cb2;
P213=Ab12'*R1*Cb1p-Cb2'*Cb1p'*Ab11'*R1*Cb1p+Q2'*Khx'*Bb1'*Cb1p-Cb2'*Cb1p'*Q1'*Khx'*Bb1'*Cb1p-Cb2'*Cb1p'*Omega2'*Khw'*Bb1'*Cb1p+Omega1'*Khw'*Bb1'*Cb1p-Cb2'*Cb1p'*R1*Ab11*Cb1p-Cb2'*Cb1p'*Bb1*Khx*Q1*Cb1p-Cb2'*Cb1p'*Bb1*Khw*Omega2*Cb1p-Cb2'*Cb1p'*R4*Cb1p;
P214=-Cb2'*Cb1p'*Bb1*Khw*Omega3-Cb2'*Cb1p'*R1*Sigma1*R';
P222=R6;
P223=Cb2'*Cb1p'*Q1'*Khx'*Bb1'*Cb1p-Q2'*Khx'*Bb1'*Cb1p-Omega1'*Khw'*Bb1'*Cb1p+Cb2'*Cb1p'*Omega2'*Khw'*Bb1'*Cb1p;
P224=zeros(n-qw,qw);
P233=Cb1p'*Ab11'*R1*Cb1p+Cb1p'*Q1'*Khx'*Bb1'*Cb1p+Cb1p'*Omega2'*Khw'*Bb1'*Cb1p+Cb1p'*R1*Ab11*Cb1p+Cb1p'*Bb1*Khx*Q1*Cb1p+Cb1p'*Bb1*Khw*Omega2*Cb1p+Cb1p'*R4*Cb1p;
P234=Cb1p'*Bb1*Khw*Omega3+Cb1p'*R1*Sigma1*R';
P244=-gama4*eye(qw,qw);

LME=list(R1-R1',R3-R3',R4-R4',R6-R6',R1*Bb1-Bb1*Rh1)
LMI=list(-([P211,P212,P213,P214;P212',P222,P223,P224;P213',P223',P233,P234;P214',P224',P234',P244]),R1,R3,R4,R6,gama4)
OBJ=[]
endfunction

R1_0=eye(qw,qw);
R3_0=eye(n-qw,n-qw);
R4_0=eye(qw,qw);
R6_0=eye(n-qw,n-qw);
Rh1_0=eye(m,m);
Khx_0=zeros(m,n);
Khw_0=zeros(m,qw);
gama4_0=1e6;

Init_guess_DRC1=list(R1_0,R3_0,R4_0,R6_0,Rh1_0,Khx_0,Khw_0,gama4_0);

Mbound=1e1;
abstol=1e-5;
nu=1;
maxiters=500;
reltol=1e-5;

//Ans_LMI=lmisolver(Init_guess_DRC1,DRC1,[Mbound,abstol,nu,maxiters,reltol]);
Ans_LMI=lmisolver(Init_guess_DRC1,DRC1);

R1=Ans_LMI(1);
R3=Ans_LMI(2);
R4=Ans_LMI(3);
R6=Ans_LMI(4);
Rh1=Ans_LMI(5);
Khx=Ans_LMI(6);
Khw=Ans_LMI(7);
gama4=Ans_LMI(8);

