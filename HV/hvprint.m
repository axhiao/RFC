clc; clear all;
% test_case='ankang20000712';
% test_case='ankang20030828';
 test_case='ankang20051001';
% test_case='ankang20100715';

str_appf = 'E:\pMOEAD for RFC1m811\HV\';
str_appf2 = strcat(str_appf,'HV_MOEAD_');
str_appf2 = strcat(str_appf2,test_case);
str_appf2 = strcat(str_appf2,'_OBJ2.txt');
HV = importdata(str_appf2);
%hv= mean (HV,1);
hv = HV(9,:);
a = 5000:5000:200000;
plot(a,hv,'-*');
set(gca,'FontSize',14,'FontName','Times New Roman');
xlabel('Evaluate times','FontSize',16,'FontName','Times New Roman');
ylabel('Hypervolume','FontSize',16,'FontName','Times New Roman');
%title('Ankang Reservoir, October 12, 2000','FontSize',16,'FontName','Times New Roman');
%title('Ankang Reservoir, August 28, 2003','FontSize',16,'FontName','Times New Roman');
 title('Ankang Reservoir, October 1st, 2005','FontSize',16,'FontName','Times New Roman');
% title('Ankang Reservoir, July 15, 2010','FontSize',16,'FontName','Times New Roman');
%set(gca,'xlim', 0:5000:200000); 

% figure;
% s = size(hv);
% b = s(2);
% boxplot(hv(:,b)');




