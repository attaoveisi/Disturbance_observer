clc
clear
close all

load SteamEng

steam = iddata([GenVolt,Speed],[Pressure,MagVolt],0.05);
steam.InputName  = {'Pressure';'MagVolt'};
steam.OutputName = {'GenVolt';'Speed'};

figure;
subplot(2,2,1)
plot(steam(:,1,1))
subplot(2,2,2)
plot(steam(:,1,2))
subplot(2,2,3)
plot(steam(:,2,1))
subplot(2,2,4)
plot(steam(:,2,2))

% A first step to get a feel for the dynamics is to look at 
% the step responses between the different channels estimated 
% directly from data:
mi = impulseest(steam,50);
clf, step(mi)

