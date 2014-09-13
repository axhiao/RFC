clc; clear all;
%str_appf = 'E:\pMOEAD for RFC2m\POS\POS_MOEAD_ankang20100715_OBJ2_RUN3.txt';
str_appf = 'E:\pMOEAD for RFC1m\UpstreamWaterLevel\UWL_MOEAD_ankang20100715_OBJ2_RUN2.txt';
POS = importdata(str_appf);
for i= 1:20
plot(POS(i,:));
hold on;
end
title('Ankang Reservoir,Shanxi,Xi''an');
xlabel('Water Level(m)');
ylabel('Reservoir Water Capacity(10^6 m^3)');