clc; clear all;
HV= importdata('E:\pMOEAD for RFC1m811\HV\HV_MOEAD_ankang20030828_OBJ2.txt');
hv= mean (HV,1);
a = 5000:5000:200000;
plot(a,hv,'-h');
xlabel('evaluate times');
ylabel('hypervolume');
title('ankang20030828');
%set(gca,'xtick', 0:5000:200000); 
%axis([0 200000])
figure;
s = size(hv);
b = s(2);
boxplot(HV(:,b)');
title('ankang20030828');










