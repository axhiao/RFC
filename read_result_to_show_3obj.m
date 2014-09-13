clc;
str_pr = 'F:\maxiaoliang\2013-1-7\¸´¼þ (2) pMOEAD\PF\';
 str_appf = 'F:\maxiaoliang\2013-1-7\¸´¼þ (2) pMOEAD\POF\';
test_case= 'UF8';
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
if strcmp(test_case, 'DTLZ1')
        point1=[[0 0 0];[0.5 0.1 0.3]];
        point2 =[[0 0 0];[0.3 0.5 0.1]];
elseif strcmp(test_case, 'DTLZ2') | strcmp(test_case, 'DTLZ3')|strcmp(test_case, 'DTLZ4')
        point1=[[0 0 0];[1 0.2 0.6]];
        point2 =[[0 0 0];[0.6 1 0.2]];
elseif strcmp(test_case, 'DTLZ6')
        point1=[[0 0 0];[0.75 0.75 5]];
        point2 =[[0 0 0];[0.75 0.13 5.3]];
elseif strcmp(test_case, 'UF8')
        point1=[[0 0 0];[1 0.2 0.6]];
        point2 =[[0 0 0];[0.6 1 0.2]];
elseif strcmp(test_case, 'UF9')
        point1=[[0 0 0];[0.65 0.1 0.4]];
        point2 =[[0 0 0];[0.1 0.65 0.4]];      
elseif strcmp(test_case, 'UF10')
        point1=[[0 0 0];[1 0.2 0.6]];
        point2 =[[0 0 0];[0.6 1 0.2]];
end
plot3(pftrue(:,1),pftrue(:,2),pftrue(:,3),'k.',Appro_pftrue(4:size(Appro_pftrue,1),1),Appro_pftrue(4:size(Appro_pftrue,1),2),Appro_pftrue(4:size(Appro_pftrue,1),3),'b.',point1(:,1),point1(:,2),point1(:,3),'g-',point2(:,1),point2(:,2),point2(:,3),'g-');
%plot3(Appro_pftrue(4:size(Appro_pftrue,1),1),Appro_pftrue(4:size(Appro_pftrue,1),2),Appro_pftrue(4:size(Appro_pftrue,1),3),'b.');
title(test_case);xlabel('f1');ylabel('f2');zlabel('f3');
end


%hold on; plot(pftrue(:,1)+0.1,pftrue(:,2)+0.1,'k-',Appro_pftrue(:,1),Appro_pftrue(:,2),'b.');

