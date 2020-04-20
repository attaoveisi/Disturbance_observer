clc
clear
close
%%
% load sys;
% A = sys.A;
% BE = sys.B;
% B = BE(:,1:2);
% E = BE(:,3);
% C = sys.C;
% DE = sys.D;
% D = DE(:,1:2);
% W = [1;0];
% 
% alfa = 0.5;
% g = 9.8;
% row = 1;
% ME = 10;
% MC = 5;
% K = 4.87;
% A = [0 1 0 0;
%      -K/ME -alfa*g K/ME 0;
%      0 0 0 1;
%      K/MC 0 -K/MC -alfa*g];
% B = [0 1/ME 0 0]';
% E = [1 -1/ME 0 0]';
% C = [1 0 0 0;
%      0 1 0 0];
% D = [0;
%      0];
% W = [0;
%      0];

load A
load B
load C
load D
load H
D=D(:,1);
load K
load L
load Kw

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
if dum4 ~= 0
%     M = rand(dum4,dum2);
    load M
    while rank([Cb1p;M]) ~= dum2;
%         M = rand(dum4,dum2);
        load M
    end
    N = [Cb1p;M];
else
    N = [Cb1p];
    M = eye(p,p);
end

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

Sai1 = Qbx*(At2-L*Ct)-A*Qbx;
Sai2 = Qbx*(Gt2+L*Ht)-A*Qby;
Sai3 = Qby;
Sai4 = Qbx*Bt2-B;
Sai5 = Qbx*Q2';

epsilon = 0.1;

[eig(At2) eig(At2-L*Ct)]
[eig(A) eig(A+B*K*0.45)]