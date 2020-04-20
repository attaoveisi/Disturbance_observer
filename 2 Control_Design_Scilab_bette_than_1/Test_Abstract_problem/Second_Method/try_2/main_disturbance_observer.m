% dumm = [];
% for kkk = 1:1000
    %% PLANT MODEL
% clc
load A
load B
load C
load H
load observer

dum1 = size(A);
n = dum1(1,1);
dum2 = size(B);
m = dum2(1,2);
dum3 = size(C);
p = dum3(1,1);
dum4 = size(H);
qw = dum4(1,2);

%% LMI Definition
setlmis([]);

delta2 = lmivar(1,[qw 1]);
beta_hat = lmivar(2,[qw n]);
gama61 = lmivar(1,[1 1]);

%LMI terms
lmiterm([1 1 1 0],1); 
lmiterm([1 1 1 beta_hat],-1,H,'s'); 
lmiterm([1 1 2 beta_hat],1,(Gama-A-Phi*C));
lmiterm([1 1 3 beta_hat],-2,(theta*C*A)); 
lmiterm([1 1 4 delta2],1,1); 
lmiterm([1 4 4 gama61],-1,1); 

lmiterm([-2 1 1 gama61],1,1);

lmiterm([-3 1 1 delta2],1,1);

LMISYS = getlmis;

%% Finding Feasible Solution 
n_dec = decnbr(LMISYS);
c = zeros(1,n_dec);
c(1:end) = 1;
options = [1e-6 1000 393 0 1];
% [copt,xopt] = mincx(LMISYS,c,options,[],[]);
[copt,xopt] = feasp(LMISYS,options);

delta2 = dec2mat(LMISYS,xopt,1);
beta_hat = dec2mat(LMISYS,xopt,2);
gama61 = dec2mat(LMISYS,xopt,3);

beta = inv(delta2)*beta_hat;

% dumm = [dumm;kkk eig(-beta*H)];
eig(-beta*H)
% kkk

save beta beta
save delta2 delta2
% end