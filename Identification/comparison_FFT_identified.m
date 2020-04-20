clc
clear
close all

%%
load Identified3_real_inputs
% load Identified

A1 = A(1:4,1:4);
B1 = B(1:4,:);
C1 = C(:,1:4);

sys = ss(A,B,C,D);
sys_reduced = ss(A1,B1,C1,D);
[sys_reduced,info1] = balancmr(sys,4);

[sv,w] = sigma(sys-sys_reduced(:,:,1));
% figure;
% loglog(w,sv,w,info1.ErrorBound(1)*ones(size(w)))
% xlabel('rad/sec');ylabel('SV');
% title('Error Bound and Model Error')

% sys=sys_reduced1;

%%
load iddata2
y=iddata2;
in=[y.Y(1).Data;y.Y(2).Data;y.Y(5).Data]';
out=[y.Y(3).Data;y.Y(4).Data]';
time = y.X.Data';
% plot(time,out)

%%
k1 = 44775;
k2 = 44775+250000;
u1=in(k1:k2,1);

k1 = 344775;
k2 = 344775+250000;
u2=in(k1:k2,2);

k1 = 644775;
k2 = 644775+250000;
u3=in(k1:k2,3);

k1 = 44775;
k2 = 44775+250000;
y11=out(k1:k2,1);
y21=out(k1:k2,2);

k1 = 344775;
k2 = 344775+250000;
y12=out(k1:k2,1);
y22=out(k1:k2,2);

k1 = 644775;
k2 = 644775+250000;
y13=out(k1:k2,1);
y23=out(k1:k2,2);
time = time(k1:k2,1);

%%
T = 1e-4;                      % Sample time
Fs = 1/T;                    % Sampling frequency
dum1 = size(time);
L = dum1(1,1);                     % Length of signal
t = (0:L-1)*T;                % Time vector
y = y11;
u = u1;

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
U = fft(u,NFFT)/L;

f = Fs/2*linspace(0,1,NFFT/2+1);
sampling_time = (1-0)/(NFFT/2+1)*(Fs/2);

figure;

subplot(2,1,1)
semilogy(f,abs(Y(1:NFFT/2+1,:)./U(1:NFFT/2+1,:)))
hold on
resp = squeeze(freqresp(sys(1,1),f,'Hz'));
resp_reduced = squeeze(freqresp(sys_reduced(1,1),f,'Hz'));
semilogy(f,abs(resp))
hold on
semilogy(f,abs(resp_reduced))
xlim([0 100])

subplot(2,1,2)
plot(f,angle((Y(1:NFFT/2+1,:)./U(1:NFFT/2+1,:)))*180/pi);
hold on
plot(f,angle(resp)*180/pi);
hold on
plot(f,angle(resp_reduced)*180/pi);
xlim([0 100])
title('TF:output 1 input1')

%%
T = 1e-4;                      % Sample time
Fs = 1/T;                    % Sampling frequency
dum1 = size(time);
L = dum1(1,1);                     % Length of signal
t = (0:L-1)*T;                % Time vector
y = y21;
u = u1;

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
U = fft(u,NFFT)/L;

f = Fs/2*linspace(0,1,NFFT/2+1);
sampling_time = (1-0)/(NFFT/2+1)*(Fs/2);

figure;

subplot(2,1,1)
semilogy(f,abs(Y(1:NFFT/2+1,:)./U(1:NFFT/2+1,:)))
hold on
resp = squeeze(freqresp(sys(2,1),f,'Hz'));
resp_reduced = squeeze(freqresp(sys_reduced(2,1),f,'Hz'));
semilogy(f,abs(resp))
hold on
semilogy(f,abs(resp_reduced))
xlim([0 100])
subplot(2,1,2)
plot(f,angle((Y(1:NFFT/2+1,:)./U(1:NFFT/2+1,:)))*180/pi);
hold on
plot(f,angle(resp)*180/pi);
hold on
plot(f,angle(resp_reduced)*180/pi);
xlim([0 100])
title('TF:output 2 input1')

%%
T = 1e-4;                      % Sample time
Fs = 1/T;                    % Sampling frequency
dum1 = size(time);
L = dum1(1,1);                     % Length of signal
t = (0:L-1)*T;                % Time vector
y = y12;
u = u2;

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
U = fft(u,NFFT)/L;

f = Fs/2*linspace(0,1,NFFT/2+1);
sampling_time = (1-0)/(NFFT/2+1)*(Fs/2);

figure;

subplot(2,1,1)
semilogy(f,abs(Y(1:NFFT/2+1,:)./U(1:NFFT/2+1,:))) 
hold on
resp = squeeze(freqresp(sys(1,2),f,'Hz'));
resp_reduced = squeeze(freqresp(sys_reduced(1,2),f,'Hz'));
semilogy(f,abs(resp))
hold on
semilogy(f,abs(resp_reduced))
xlim([0 100])
subplot(2,1,2)
plot(f,angle((Y(1:NFFT/2+1,:)./U(1:NFFT/2+1,:)))*180/pi);
hold on
plot(f,angle(resp)*180/pi);
hold on
plot(f,angle(resp_reduced)*180/pi);
xlim([0 100])
title('TF:output 1 input2')

%%
T = 1e-4;                      % Sample time
Fs = 1/T;                    % Sampling frequency
dum1 = size(time);
L = dum1(1,1);                     % Length of signal
t = (0:L-1)*T;                % Time vector
y = y22;
u = u2;

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
U = fft(u,NFFT)/L;

f = Fs/2*linspace(0,1,NFFT/2+1);
sampling_time = (1-0)/(NFFT/2+1)*(Fs/2);

figure;

subplot(2,1,1)
semilogy(f,abs(Y(1:NFFT/2+1,:)./U(1:NFFT/2+1,:))) 
hold on
resp = squeeze(freqresp(sys(2,2),f,'Hz'));
resp_reduced = squeeze(freqresp(sys_reduced(2,2),f,'Hz'));
semilogy(f,abs(resp))
hold on
semilogy(f,abs(resp_reduced))
xlim([0 100])
subplot(2,1,2)
plot(f,angle((Y(1:NFFT/2+1,:)./U(1:NFFT/2+1,:)))*180/pi);
hold on
plot(f,angle(resp)*180/pi);
hold on
plot(f,angle(resp_reduced)*180/pi);
xlim([0 100])
title('TF:output 2 input2')

%%
T = 1e-4;                      % Sample time
Fs = 1/T;                    % Sampling frequency
dum1 = size(time);
L = dum1(1,1);                     % Length of signal
t = (0:L-1)*T;                % Time vector
y = y13;
u = u3;

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
U = fft(u,NFFT)/L;

f = Fs/2*linspace(0,1,NFFT/2+1);
sampling_time = (1-0)/(NFFT/2+1)*(Fs/2);

figure;

subplot(2,1,1)
semilogy(f,abs(Y(1:NFFT/2+1,:)./U(1:NFFT/2+1,:))) 
hold on
resp = squeeze(freqresp(sys(1,3),f,'Hz'));
resp_reduced = squeeze(freqresp(sys_reduced(1,3),f,'Hz'));
semilogy(f,abs(resp))
hold on
semilogy(f,abs(resp_reduced))
xlim([0 100])
subplot(2,1,2)
plot(f,angle((Y(1:NFFT/2+1,:)./U(1:NFFT/2+1,:)))*180/pi);
hold on
plot(f,angle(resp)*180/pi);
hold on
plot(f,angle(resp_reduced)*180/pi);
xlim([0 100])
title('TF:output 1 input3')

%%
T = 1e-4;                      % Sample time
Fs = 1/T;                    % Sampling frequency
dum1 = size(time);
L = dum1(1,1);                     % Length of signal
t = (0:L-1)*T;                % Time vector
y = y23;
u = u3;

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
U = fft(u,NFFT)/L;

f = Fs/2*linspace(0,1,NFFT/2+1);
sampling_time = (1-0)/(NFFT/2+1)*(Fs/2);

figure;

subplot(2,1,1)
semilogy(f,abs(Y(1:NFFT/2+1,:)./U(1:NFFT/2+1,:))) 
hold on
resp = squeeze(freqresp(sys(2,3),f,'Hz'));
resp_reduced = squeeze(freqresp(sys_reduced(2,3),f,'Hz'));
semilogy(f,abs(resp))
hold on
semilogy(f,abs(resp_reduced))
xlim([0 100])
subplot(2,1,2)
plot(f,angle((Y(1:NFFT/2+1,:)./U(1:NFFT/2+1,:)))*180/pi);
hold on
plot(f,angle(resp)*180/pi);
hold on
plot(f,angle(resp_reduced)*180/pi);
xlim([0 100])
title('TF:output 2 input3')