clear 
close all
clc

% m1=15;
% m2=1+7.8;
% m3=43.4;
% k1=31000;
% k2=18000;
% k3=44130;
% c1=830;
% c2=200;
% c3=1485;
% 
% A=[0 1 0 0 0 0;
%     -k1/m1 -(c1+c2)/m1 k2/m1 c2/m1 0 0;
%     0 -1 0 1 0 0;
%     0 c2/m2 -k2/m2 -(c2+c3)/m2 k3/m2 c3/m2;
%     0 0 0 -1 0 1;
%     0 0 0 c3/m3 -k3/m3 -c3/m3];
% 
% B=eye(6,3);
% 
% H=[-1;c1/m1;0;0;0;0];
% 
% C1n=[0 0 0 c3/m3 -k3/m3 -c3/m3];  %%% Output for Ride Comfort
% 
% C2n=[1 0 0 0 0 0];  %%% Output for Suspension Deflection
% 
% C=[C2n;C1n;0 0 1 0 0 0];

nout = 3;
nin = 3;
nstate = 6;
for k = 1:100000
    ii = 1;
    sys_rand = rss(nstate,nout,nin);

    A = sys_rand.A;
    B = sys_rand.B;
    H = B(:,end);
    B = B(:,1:end-1);
    C = sys_rand.C;

    sys_dum = tf(sys_rand);
    dumm = zeros(nstate,nin*nout);
    for i = 1:nin
        for j = 1:nout
            dummm = squeeze(zero(sys_dum(j,i)));
            [x1,x2] = size(dummm);
            if x1 < nstate
                dummm = [dummm;-ones(nstate-x1,1)];
            end
            dumm(:,ii) = dummm;
            ii = ii+1;
        end
    end
    if real(dumm) < 0
        real(pole(sys_rand))
        real(dumm)
        break
    end
end
        
%%
save A A
save B B
save H H 
save C C

%% PLANT MODEL
load A
load B
load C
load H

dum1 = size(A);
n = dum1(1,1);
dum2 = size(B);
m = dum2(1,2);
dum3 = size(C);
p = dum3(1,1);
dum4 = size(H);
qw = dum4(1,2);

CHp=(C*H)'*((C*H)'*(C*H))^-1;

theta1 = -H*CHp;
theta2 = eye(p,p)-(C*H)*CHp;
eta1 = eye(n,n)-H*CHp*C;
eta2 = (eye(p,p)-(C*H)*CHp)*C;

%% Finding Y through GA optimization
nvars = n*p;
lb = -1e3*ones(1,n);
ub = 1e3*ones(1,n);
PopulationSize_Data = 100;
EliteCount_Data = 2;
CrossoverFraction_Data = 0.9;
MigrationFraction_Data = 0.25;
Generations_Data = 2000;
StallGenLimit_Data = 100;
TolFun_Data = 1e-10;
TolCon_Data = 1e-10;
[Y_dum,fval,exitflag,output,population,score] = optimization_for_Y(nvars,lb,ub,PopulationSize_Data,EliteCount_Data,CrossoverFraction_Data,MigrationFraction_Data,Generations_Data,StallGenLimit_Data,TolFun_Data,TolCon_Data);

Y = [Y_dum(1) Y_dum(2) Y_dum(3);
     Y_dum(4) Y_dum(5) Y_dum(6);
     Y_dum(7) Y_dum(8) Y_dum(9);
     Y_dum(10) Y_dum(11) Y_dum(12);
     Y_dum(13) Y_dum(14) Y_dum(15);
     Y_dum(16) Y_dum(17) Y_dum(18)];

eta = eta1+Y*eta2;
theta = theta1+Y*theta2;
Delta = (eye(n,n)+(-H*CHp+Y*(eye(p,p)-(C*H)*CHp))*C)*B;

eta*H
close all

%% LMI Definition
setlmis([]);

X_hat = lmivar(2,[n p]);
Phi_hat = lmivar(2,[n p]);
delta1 = lmivar(1,[n 1]);
% gama5 = lmivar(1,[1 1]);

%LMI terms
lmiterm([1 1 1 delta1],A'*eta',1,'s'); 
lmiterm([1 1 1 X_hat],-1,C,'s'); 
lmiterm([1 1 1 Phi_hat],-1,C,'s'); 
% lmiterm([1 2 2 0],C'*C);
% lmiterm([1 3 3 gama5],-1,1); 

lmiterm([-2 1 1 delta1],1,1); 
lmiterm([-2 1 1 0],-0.01); 

% lmiterm([-3 1 1 gama5],1,1);

LMISYS = getlmis;

%% Finding Feasible Solution 
n_dec = decnbr(LMISYS);
c = zeros(1,n_dec);
c(1:end) = 1;
options = [1e-6 1000 0 0 0];
% [copt,xopt] = mincx(LMISYS,c,options,[],[]);
[copt,xopt] = feasp(LMISYS);

X_hat = dec2mat(LMISYS,xopt,1);
Phi_hat = dec2mat(LMISYS,xopt,2);
delta1 = dec2mat(LMISYS,xopt,3);
% gama5 = dec2mat(LMISYS,xopt,4);

%% Observer parameters
X = inv(delta1)*X_hat;
Phi = inv(delta1)*Phi_hat;
Gama = eta*A-X*C;

eig(Gama-Phi*C)
eig(Gama)

save('observer','Gama','Phi','Delta','theta');
SimOut = sim('main_observer_sim');