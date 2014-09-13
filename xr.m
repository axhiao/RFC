 %16\14\18\24
 clc;clear all;
 a=6:6:6*24;
str_appf = 'E:\421\pMOEAD for RFC\515\PUpstreamWaterLevel\';
%test_case='ankang20000712';
%test_case='ankang20030828';
%test_case='ankang20051001';
test_case='ankang20100715';
count = 1;

for i=1:count
str_appf2 = strcat(str_appf, 'UWL_MOEAD_');
str_appf2 = strcat(str_appf2, test_case);
str_appf2 = strcat(str_appf2, '_OBJ2_RUN');
str_appf2 = strcat(str_appf2, num2str(i));
str_appf2 = strcat(str_appf2, '.txt');
PUML = importdata(str_appf2);
%figure(10+floor(rand*1000));

for j=1:20
plot(a,PUML(j,:),'-');
h(j)=max(PUML(j,:));
hold on;
%figure
%plot(j,h(j),'*')
end
title(test_case);xlabel('Time');ylabel('PUpstreamwaterlevel');
end
figure
b=1:1:20;
plot(b,h,'o') 

hold on;
%line([0,20],[325,325],'r-');
plot([0 20], [325 325],'r-')
grid on;
title(test_case);xlabel('number');ylabel('maxPUpstreamwaterlevel');
figure
 
 
 str_appf = 'E:\421\pMOEAD for RFC\515\PUpstreamWaterLevel\';
%test_case='ankang20000712';
%test_case='ankang20030828';
%test_case='ankang20051001';
test_case='ankang20100715';
count = 1;

for i=1:count
str_appf2 = strcat(str_appf, 'UWL_MOEAD_');
str_appf2 = strcat(str_appf2, test_case);
str_appf2 = strcat(str_appf2, '_OBJ2_RUN');
str_appf2 = strcat(str_appf2, num2str(i));
str_appf2 = strcat(str_appf2, '.txt');
UML = importdata(str_appf2);
%figure(10+floor(rand*1000));

for j=1:20
plot(a,UML(j,:));
g(j)=max(UML(j,:));
hold on;
%figure
%plot(j,h(j),'*')
end
title(test_case);xlabel('Time');ylabel('Upstreamwaterlevel');
end
figure
b=1:1:20;
plot(b,h,'bO',b,g,'mH') 
legend('newpopulation','population')
%legend('boxoff');

set(gca,'xtick',[0:1:20]);   %…Ë÷√x÷·
set(gca,'ytick',[315:1:332]);   %…Ë÷√y÷·   
hold on;
title(test_case);xlabel('number');ylabel('maxUpstreamwaterlevel');
%line([0,20],[325,325],'r-');
 plot([0 20], [325 325],'r-')
 grid on;
 