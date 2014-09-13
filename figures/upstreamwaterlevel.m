clc; clear all;
str_appf = 'E:\pMOEAD for RFC1m\UpstreamWaterLevel\UWL_MOEAD_ankang20001012_OBJ2_RUN3.txt';
UWL = importdata(str_appf);
for i= 1:20
plot(UWL(i,:));
hold on;
end
title('Ankang Reservoir, October 12, 2000');
set(gca,'xtick', 0:5:20); %???????????????????????????????????????????????
%6 3 4 6
xlabel('Dispatching Periods (Time Interval: 6 Hours)');
ylabel('Upstream Water Level (m)')

% clc; clear all;
% str_appf = 'E:\pMOEAD for RFC1m\UpstreamWaterLevel\UWL_MOEAD_ankang20030828_OBJ2_RUN3.txt';
% UWL = importdata(str_appf);
% for i= 1:20
% plot(UWL(i,:));
% hold on;
% end
% title('Ankang Reservoir, August 28, 2003');
% set(gca,'xtick', 0:5:15); 
% %6 3 4 6
% xlabel('Dispatching Periods (Time Interval: 3 Hours)');
% ylabel('Upstream Water Level (m)')

% clc; clear all;
% str_appf = 'E:\pMOEAD for RFC1m\UpstreamWaterLevel\UWL_MOEAD_ankang20051001_OBJ2_RUN3.txt';
% UWL = importdata(str_appf);
% for i= 1:20
% plot(UWL(i,:));
% hold on;
% end
% title('Ankang Reservoir, October 1st, 2005');
% set(gca,'xtick', 0:5:25); 
% %6 3 4 6
% xlabel('Dispatching Periods (Time Interval: 4 Hours)');
% ylabel('Upstream Water Level (m)')

% clc; clear all;
% str_appf = 'E:\pMOEAD for RFC1m\UpstreamWaterLevel\UWL_MOEAD_ankang20100715_OBJ2_RUN3.txt';
% UWL = importdata(str_appf);
% for i= 1:20
% plot(UWL(i,:));
% hold on;
% end
% title('Ankang Reservoir, July 15, 2010');
% set(gca,'xtick', 0:5:25); 
% %6 3 4 6
% xlabel('Dispatching Periods (Time Interval: 6 Hours)');
% ylabel('Upstream Water Level (m)')