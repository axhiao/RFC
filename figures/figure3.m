clc; clear all;
str_appf = 'E:\pMOEAD for RFCfu\Qdata\ankang20001012.txt';
Inflow = importdata(str_appf);
plot(Inflow(:,1));
hold on;
plot(Inflow(:,1),'.');
title('Ankang Reservoir, October 12, 2000');
xlabel('Periods Index(Time Interval: 1 Hour)');
ylabel('Reservoir Inflow Volume(m^3/s)')

% clc; clear all;
% str_appf = 'E:\pMOEAD for RFCfu\Qdata\ankang20030828.txt';
% Inflow = importdata(str_appf);
% plot(Inflow(:,1));
% hold on;
% plot(Inflow(:,1),'.');
% title('Ankang Reservoir, August 28, 2003');
% xlabel('Periods Index(Time Interval: 1 Hour)');
% ylabel('Reservoir Inflow Volume(m^3/s)')
% axis([0 50 0 14000]);

% clc; clear all;
% str_appf = 'E:\pMOEAD for RFCfu\Qdata\ankang20051001.txt';
% Inflow = importdata(str_appf);
% plot(Inflow(:,1));
% hold on;
% plot(Inflow(:,1),'.');
% title('Ankang Reservoir, October 1, 2005');
% xlabel('Periods Index(Time Interval: 1 Hour)');
% ylabel('Reservoir Inflow Volume(m^3/s)')

% clc; clear all;
% str_appf = 'E:\pMOEAD for RFCfu\Qdata\ankang20100715.txt';
% Inflow = importdata(str_appf);
% plot(Inflow(:,1));
% hold on;
% plot(Inflow(:,1),'.');
% title('Ankang Reservoir,July 15, 2010');
% xlabel('Periods Index(Time Interval: 1 Hour)');
% ylabel('Reservoir Inflow Volume(m^3/s)')