str_appf = 'perfectPareto\';
%test_case='ankang20000712';
test_case='ankang20030828';
%test_case='ankang20051001';
%test_case='ankang20100715';
str_appf2 = strcat(str_appf, 'PF_');
str_appf2 = strcat(str_appf2, test_case);
str_appf2 = strcat(str_appf2, '.txt');
PF = importdata(str_appf2);
plot(PF(:,1),PF(:,2),'.')
hold on


str_appf = 'POF\';
%str_appf = 'PPOF\';
%test_case='ankang20000712';
%test_case='ankang20030828';
%test_case='ankang20051001';
%test_case='ankang20100715';
count = 10;

%POF_MOEAD_ankang20000712_OBJ2_RUN1
for i=10:count
str_appf2 = strcat(str_appf, 'POF_MOEAD_');
str_appf2 = strcat(str_appf2, test_case);
str_appf2 = strcat(str_appf2, '_OBJ2_RUN');
str_appf2 = strcat(str_appf2, num2str(i));
str_appf2 = strcat(str_appf2, '.txt');
POF = importdata(str_appf2);
%figure(10+floor(rand*1000));
%set(gca,'xtick',200:5:335);
%set(gca,'ytick',5000:500:10000)
%axis ([200,335,4000,10000])
plot(POF(:,1),POF(:,2),'rH')
hold on;
end
%axis ([200,400,6000,20000])
title(test_case);xlabel('Maximum Upstream Water Level (m)');ylabel('Maximum Out Flow (m3/s)');



count=10;
str_appf = 'PPOF\';
for i=10:count
str_appf2 = strcat(str_appf, 'POF_MOEAD_');
str_appf2 = strcat(str_appf2, test_case);
str_appf2 = strcat(str_appf2, '_OBJ2_RUN');
str_appf2 = strcat(str_appf2, num2str(i));
str_appf2 = strcat(str_appf2, '.txt');
PPOF = importdata(str_appf2);
%figure(10+floor(rand*1000));
%set(gca,'xtick',200:5:335);
%set(gca,'ytick',5000:500:10000)
%axis ([200,335,4000,10000])
plot(PPOF(:,1),PPOF(:,2),'bO')
hold on;
end
