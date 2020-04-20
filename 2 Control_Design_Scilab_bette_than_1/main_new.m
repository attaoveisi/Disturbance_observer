clc
clear
close

load A
load B
load C
load H
load M
load N
load J
load L

Ap = A;
Bp = B;
Cp = C;
Hp = H;

m = rank(B);
n = rank(A);
qw = rank(H);
p = rank(C);

%%
[Q,Sigmaw,R] = svd(H);
Sigma1 = Sigmaw(1:qw,:);
Q1 = Q(:,1:qw);
Q2 = Q(:,qw+1:end);
%%
Ab = Q'*A*Q;
Bb = Q'*B;
Ab11 = Ab(1:qw,1:qw);
Ab12 = Ab(1:qw,qw+1:end);
Ab21 = Ab(qw+1:end,1:qw);
Ab22 = Ab(qw+1:end,qw+1:end);
Bb1 = Bb(1:qw,:);
Bb2 = Bb(qw+1:end,:);
Cb = C*Q;
Cb1 = Cb(:,1:qw);
Cb2 = Cb(:,qw+1:end);
Cb1p = (Cb1'*Cb1)^-1*Cb1';
dum1 = size(Cb1p);
dum2 = dum1(1,2);
dum3 = dum1(1,1);
dum4 = dum2-dum3;

At2 = Ab22-Ab21*Cb1p*Cb2;
Bt2 = Bb2;
Gt2 = Ab21*Cb1p;
Ct = M*(eye(p)-Cb1*Cb1p)*Cb2;
Ht = M*(eye(p)-Cb1*Cb1p);
dum5 = size(Cb1p*Cb2);
dum6 = dum5(1,2);
Qbx = Q*[-Cb1p*Cb2;eye(dum6)];
dum7 = size(Cb1p);
dum8 = dum7(1,2);
dum9 = size(Q);
dum10 = dum9(1,2);
dum11 = dum10-dum7(1,1);
Qby = Q*[Cb1p;zeros(dum11,dum8)];
dum15=size(M);

eig(N-L*Ct)
