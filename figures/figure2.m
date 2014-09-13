clc; clear all;
str_appf = 'E:\pMOEAD for RFC1m\Qdata\L2C.txt';
L2C = importdata(str_appf);
plot(L2C(:,1),L2C(:,2));
title('Ankang Reservoir,Shanxi,Xi''an');
xlabel('Water Level(m)');
ylabel('Reservoir Water Capacity(10^6 m^3)');
