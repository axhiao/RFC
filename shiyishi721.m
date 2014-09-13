clc; clear all;
PF = importdata('E:\pMOEAD for RFC1m\perfectPareto\PF_ankang20001012.txt');
plot(PF(:,1),PF(:,2),'b.')
hold on
POF = importdata('E:\pMOEAD for RFC1m\POF\POF_MOEAD_ankang20000712_OBJ2_RUN5.txt');
plot(POF(:,1),POF(:,2),'rO')


% str_appf = 'E:\pMOEAD for RFCfu\Qdata\L2C.txt';
% L2C = importdata(str_appf);
% plot(L2C(:,1),L2C(:,2));
% title('Ankang Reservoir,Shanxi,Xi''an');
% xlabel('Water Level(m)');
% ylabel('Reservoir Water Capacity(10^6 m^3)');
