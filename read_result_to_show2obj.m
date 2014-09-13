clc;
str_pr = 'F:\maxiaoliang\2013-1-7\¸´¼þ (2) pMOEAD\PF\';
 str_appf = 'F:\maxiaoliang\2013-1-7\¸´¼þ (2) pMOEAD\POF\';
test_case= 'UF7';
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
if strcmp(test_case, 'ZDT1') | strcmp(test_case, 'ZDT2')  | strcmp(test_case, 'ZDT3') ...
        | strcmp(test_case, 'ZDT4') | strcmp(test_case, 'ZDT6') |strcmp(test_case, 'UF1') ...
        | strcmp(test_case, 'UF2')  |strcmp(test_case, 'UF3')  | strcmp(test_case, 'UF4') ...
        | strcmp(test_case, 'UF5') | strcmp(test_case, 'UF7')
    point=[[0 0];[1 0.5]];
    point2=[[0 0];[0.5 1]];
elseif strcmp(test_case, 'UF6')
    point=[[0 0];[0.6 1]];
    point2=[[0 0];[1 0.1428]];
end
plot(pftrue(:,1),pftrue(:,2),'k-',Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.',point(:,1),point(:,2),'g-',point2(:,1),point2(:,2),'g-');
title(test_case);xlabel('f1');ylabel('f2');
end


%hold on; plot(pftrue(:,1)+0.1,pftrue(:,2)+0.1,'k-',Appro_pftrue(:,1),Appro_pftrue(:,2),'b.');

