%%%%%%% mu Synthesis

clear all;
clc;
%%% System Parameters
m1=15;
m2=1+7.8;
m3=43.4;
k1=31000;
k2=18000;
k3=44130;
c1=830;
c2=200;
c3=1485;
%% System Dynamic Model by Considering Uncertainties
delta=.02;

a11=0;
a12=1;
a13=0;
a14=0;
a15=0;
a16=0;
a21=-k1/m1;
a22=-(c1+c2)/m1;
a23=k2/m1;
a24=c2/m1;
a25=0;
a26=0;
a31=0;
a32=-1;
a33=0;
a34=1;
a35=0;
a36=0;
a41=0;
a42=ureal('a42',c2/m2,'PlusMinus',[-delta^2 delta^2]);
a43=ureal('a43',-k2/m2,'PlusMinus',[-delta^2 delta^2]);
a44=ureal('a44',-(c2+c3)/m2,'PlusMinus',[-delta^2 delta^2]);
a45=ureal('a45',k3/m2,'PlusMinus',[-delta^2 delta^2]);
a46=ureal('a46',c3/m2,'PlusMinus',[-delta^2 delta^2]);
a51=0;
a52=0;
a53=0;
a54=-1;
a55=0;
a56=1;
a61=0;
a62=0;
a63=0;
a64=ureal('a64',c3/m3,'PlusMinus',[-delta^2 delta^2]);
a65=ureal('a65',-k3/m3,'PlusMinus',[-delta^2 delta^2]);
a66=ureal('a66',-c3/m3,'PlusMinus',[-delta^2 delta^2]);

A=[a11 a12 a13 a14 a15 a16;
    a21 a22 a23 a24 a25 a26;
    a31 a32 a33 a34 a35 a36;
    a41 a42 a43 a44 a45 a46;
    a51 a52 a53 a54 a55 a56;
    a61 a62 a63 a64 a65 a66;];

b1=0;
b2=-1/m1;
b3=0;
b4=ureal('b4',0,'PlusMinus',[-delta^2 delta^2]);
b5=0;
b6=0;

B=[b1;b2;b3;b4;b5;b6];

c11=0;
c12=0;
c13=0;
c14=ureal('c14',c3/m3,'PlusMinus',[-delta^2 delta^2]);
c15=ureal('c15',-k3/m3,'PlusMinus',[-delta^2 delta^2]);
c16=ureal('c16',-c3/m3,'PlusMinus',[-delta^2 delta^2]);

Bw=[-1;c1/m1;0;0;0;0];

C1=[c11 c12 c13 c14 c15 c16];  %%% Output for Ride Comfort

C2=[1 0 0 0 0 0];  %%% Output for Suspension Deflection

cy11=ureal('cy11',1,'PlusMinus',[-delta^2 delta^2]);
cy23=ureal('cy23',1,'PlusMinus',[-delta^2 delta^2]);

Cy=[cy11 0 0 0 0 0;0 0 cy23 0 0 0];   %%% Output for Feedback

%% Nominal System
An=[0 1 0 0 0 0;
    -k1/m1 -(c1+c2)/m1 k2/m1 c2/m1 0 0;
    0 -1 0 1 0 0;
    0 c2/m2 -k2/m2 -(c2+c3)/m2 k3/m2 c3/m2;
    0 0 0 -1 0 1;
    0 0 0 c3/m3 -k3/m3 -c3/m3];

Bn=[0;-1/m1;0;0;0;0];

Bwn=[-1;c1/m1;0;0;0;0];

C1n=[0 0 0 c3/m3 -k3/m3 -c3/m3];  %%% Output for Ride Comfort

C2n=[1 0 0 0 0 0];  %%% Output for Suspension Deflection

Cyn=[1 0 0 0 0 0;0 0 1 0 0 0];%;0 0 0 0 1 0];

sysnom=ss(An,[Bwn Bn],[C1n;C2n;Cyn],0);

%% Delay
tau=ureal('tau',7.5e-3,'PlusMinus',[-7.5e-3 7.5e-3]);

w_tau=tf([-tau/2 1],[tau/2 1]);

%%%%%%% Disturbance Formulation

a=.1;    % meter
l=2;     % meter
V0=30*(1000/3600);   % m/s

t1=0:.01:l/V0;
t2=l/V0+.01:.01:3;
t=[t1 t2];

twn=0:.01:1;

z01=a/2*(1-cos(2*pi*V0*t1./l));
z02=0*t2;
z0=[z01 z02];

w1=a*pi*V0/l*sin(2*pi*V0*t1./l);
w2=0*t2;
w=[w1 w2];

%%%%%% Augmented Plant for mu
sys1=ss(A,[Bw B],[C1;C2;Cy],0);

systemnames = 'sys1 w_tau'; 
inputvar = '[dist;fs]'; 
outputvar = '[750*sys1(1);fs;sys1(3);sys1(4)]'; 
input_to_sys1 = '[dist;w_tau]'; 
input_to_w_tau = '[fs]'; 
sys2 = sysic; 

%%%%%%%%  D-K syn
nmeas=2;
ncont=1;

[Kmu,clp,bnd] = dksyn(sys2,nmeas,ncont);

%%%%%%%%%%%%%%%% Augmented Plant for Hinf
% systemnames = 'sysnom w_tau'; 
% inputvar = '[dist;fs]'; 
% outputvar = '[100*sysnom(1);fs;sysnom(3);sysnom(4)]'; 
% input_to_sysnom = '[dist;w_tau]'; 
% input_to_w_tau = '[fs]'; 
% sys2 = sysic; 
% 
% %%%%%%%% Hinf syn
% nmeas=2;
% ncont=1;
% 
% [Khinf,CL,GAM,INFO] = hinfsyn(sys2,nmeas,ncont);

%%%%% Closed Loop for Hinf
%D_new=[0 0;0 0;0 1;0 0;0 0];
%sys4hinf=ss(An,[Bwn Bn],[C1n;C2n;0 0 0 0 0 0;Cyn],D_new);

%CLhinf=lft(sys4hinf,Khinf);

%%%%% Closed Loop for Mu
D_new=[0 0;0 0;0 1;0 0;0 0];
sys4mu=ss(An,[Bwn Bn],[C1n;C2n;0 0 0 0 0 0;Cyn],D_new);

CLmu=lft(sys4mu,Kmu);

%%%%%%%%%% Simulation

%lsim(usample(CLmu,10),w,t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Uncertain System
sysuncert=ss(A,[Bw B],[C1;Cy],0);

systemnames = 'sysuncert w_tau Kmu'; 
inputvar = '[dist]'; 
outputvar = '[sysuncert(1)]'; 
input_to_sysuncert = '[dist;w_tau]'; 
input_to_w_tau = '[Kmu]'; 
input_to_Kmu = '[sysuncert(2);sysuncert(3)]'; 
CLuncert = sysic; 


[perfmarg,perfmargunc,report,info] = robustperf(CLuncert);
semilogx(info.MussvBnds)
xlabel('Frequency (rad/sec)');
ylabel('Mu upper/lower bounds');