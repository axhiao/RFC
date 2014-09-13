str_pr = 'PF\';
 str_appf = 'POF\';
test_case= 'ZDT1';
count = 3;

str_pr = strcat(str_pr, test_case);
str_pr = strcat(str_pr, '.dat');
pftrue = importdata(str_pr);

for i=1:count
str_appf2 = strcat(str_appf, 'POF_MOEAD_');
str_appf2 = strcat(str_appf2, test_case);
str_appf2 = strcat(str_appf2, '_RUN');
str_appf2 = strcat(str_appf2, num2str(i));
str_appf2 = strcat(str_appf2, '.txt');
Appro_pftrue = importdata(str_appf2);

figure(10+floor(rand*1000));
point=[[0 0];[1 1]];
plot(pftrue(:,1),pftrue(:,2),'k.',Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.',point(:,1),point(:,2),'g-');
%plot(Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.')
title(test_case);xlabel('f1');ylabel('f2');
end