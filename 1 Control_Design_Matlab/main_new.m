clc
clear
close

load A
load B
load C
load D
load H

Ap = A;
Bp = B;
Cp = C;
Dp = D;
Hp = H;

global A1

A=A(1:4,1:4);
B=B(1:4,:);
C=C(:,1:4);
H=H(1:4,:);
A1 = A;

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
    M = rand(dum4,dum2);
%     load M
    while rank([Cb1p;M]) ~= dum2;
        M = rand(dum4,dum2);
%         load M
    end
    N = [Cb1p;M];
else
    N = Cb1p;
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
dum15=size(M);

%%
global Kh
setlmis([]);

P1 = lmivar(1,[n-qw 1]);
P2 = lmivar(1,[n 1]);
Kh = lmivar(2,[m n]);
Lh = lmivar(2,[n-qw dum15(1,1)]);
gama = lmivar(1,[1 0]);

%LMI terms
lmiterm([1 1 1 P2],A',1,'s'); % LMI #1: A'P2+P2A
lmiterm([1 1 1 Kh],B,1,'S'); % LMI #1: B*Kh+Kh'*B'
lmiterm([1 1 1 0],C'*C); % LMI #1:C'*C
lmiterm([1 1 2 Kh],-B,1); % LMI #1: -B*Kh
lmiterm([1 1 3 P2],1,H); % LMI #1: P2H
lmiterm([1 1 4 P2],1,B); % LMI #1: P2*B
lmiterm([1 2 2 P1],Q2*At2',Q2','s'); % LMI #1: 
lmiterm([1 2 2 Lh],-Q2,Ct*Q2','s'); % LMI #1: 
lmiterm([1 2 2 0],Q2*Q2'); % LMI #1: 
lmiterm([1 3 3 0],-gama); % LMI #1: -gama*I
lmiterm([1 4 4 0],0); % LMI #1: 0

lmiterm([-2 1 1 P1],1,1); % LMI #2: P>0

lmiterm([-3 1 1 P2],1,1); % LMI #3: R>0

lmiterm([-4 1 1 gama],1,1); % LMI #4: eps1>0

% lmiterm([-9 1 1 gama_s],1,1); % LMI #9: gama_s>0
% % lmiterm([-9 1 1 0],-0.95); % LMI #9: gama_s>0.95
% 
% lmiterm([10 1 1 gama_s],1,1); % LMI #9: gama_s<1
% lmiterm([10 1 1 0],-1); % LMI #9: gama_s<1

LMISYS = getlmis;

%% Finding Feasible Solution 
ind_c = decnbr(LMISYS);
opt_c = zeros(1,ind_c);
opt_c(end) = 1;
% [copt,xopt] = mincx(LMISYS,opt_c,[0 1000 -1 99 0]);
[~,xopt] = feasp(LMISYS,[0 1000 -1 10 0]);
% [~,xopt] = feasp(LMISYS);
P1 = dec2mat(LMISYS,xopt,1);
P2 = dec2mat(LMISYS,xopt,2);
Kh = dec2mat(LMISYS,xopt,3);
Lh = dec2mat(LMISYS,xopt,4);
gama = dec2mat(LMISYS,xopt,5);

save('P2','P2');

%% PB-BPh
global PP
global BB
PP = P2;
BB = B;
nvars = m^2;
lb = -1*ones(1,nvars);
ub = 1*ones(1,nvars);
pop_size = 1000;
Gen = 2000;
options = gaoptimset; 
options = gaoptimset(options,'PopulationType','doubleVector');
options = gaoptimset(options,'PopInitRange',[lb;ub]);
options = gaoptimset(options,'PopulationSize',pop_size);
options = gaoptimset(options,'EliteCount',2);
options = gaoptimset(options,'CrossoverFraction',0.6);
options = gaoptimset(options,'MigrationDirection','forward');
options = gaoptimset(options,'MigrationInterval',15);
options = gaoptimset(options,'MigrationFraction',0.1);
options = gaoptimset(options,'Generations',Gen);
options = gaoptimset(options,'TimeLimit',Inf);
options = gaoptimset(options,'FitnessLimit',-Inf);
options = gaoptimset(options,'StallGenLimit',660);
options = gaoptimset(options,'StallTimeLimit',Inf);
options = gaoptimset(options,'TolFun',1.0000e-06);
options = gaoptimset(options,'TolCon',1.0000e-06);
options = gaoptimset(options,'InitialPopulation',[]);
%options = gaoptimset(options,'InitialPopulation',x_init);
options = gaoptimset(options,'InitialScores',[]);
options = gaoptimset(options,'InitialPenalty',10);
options = gaoptimset(options,'PenaltyFactor',100);
options = gaoptimset(options,'PlotInterval',1);
options = gaoptimset(options,'CreationFcn',@gacreationuniform);
options = gaoptimset(options,'FitnessScalingFcn',@fitscalingprop);
options = gaoptimset(options,'SelectionFcn',@selectionroulette);
options = gaoptimset(options,'CrossoverFcn',@crossoverheuristic);
%options = gaoptimset(options,'MutationFcn', @mutationgaussian);
% options = gaoptimset(options,'MutationFcn', {@mutationuniform,.05});
options = gaoptimset(options,'MutationFcn', {@mutationadaptfeasible,.1});
options = gaoptimset(options,'Display','off');
options = gaoptimset(options,'PlotFcns',@gaplotbestf);
options = gaoptimset(options,'OutputFcns',[]);
options = gaoptimset(options,'Vectorized','off') ;
options = gaoptimset(options,'FitnessScalingFcn',@fitscalingrank);

[x_sol,FVAL,REASON,OUTPUT,POPULATION,SCORES] = ga(@PhFunc_ga,nvars,[],[],[],[],lb,ub,[],options);

% 4 inputs

% Ph = [x_sol(1,1) x_sol(1,2) x_sol(1,3) x_sol(1,4);...
%     x_sol(1,2) x_sol(1,5) x_sol(1,6) x_sol(1,7);...
%     x_sol(1,3) x_sol(1,6) x_sol(1,8) x_sol(1,9);
%     x_sol(1,4) x_sol(1,7) x_sol(1,9) x_sol(1,10)];

% 2 inputs

Ph2 = [x_sol(1,1) x_sol(1,2);...
    x_sol(1,3) x_sol(1,4)];

K=Ph2^-1*Kh;
L=P1^-1*Lh;

Sai1 = Qbx*(At2-L*Ct)-A*Qbx;
Sai2 = Qbx*(Gt2+L*Ht)-A*Qby;
Sai3 = Qby;
Sai4 = Qbx*Bt2-B;
Sai5 = Qbx*Q2';


% [eig(At2) eig(At2-L*Ct)];