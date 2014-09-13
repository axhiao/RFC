clc;clear all;close all;
%16\14\18\24
%T=6,3,4,6
clc;clear all;
a=6:6:6*16;
test_case='ankang20000712';

% a=3:3:3*14;
% test_case='ankang20030828';

% a=4:4:4*18;
% test_case='ankang20051001';

% a=6:6:6*24;
% test_case='ankang20100715';

str_appf = 'E:\pMOEAD for RFCfu\UpstreamWaterLevelmid\';

%30000_UWL_MOEAD_ankang20100715_OBJ2_RUN1

%count = 10;
%for i=21000:1500:49500
for i=84000:6000:198000
str_appf2 = strcat(str_appf, num2str(i));
str_appf2 = strcat(str_appf2, '_UWL_MOEAD_');
str_appf2 = strcat(str_appf2, test_case);
str_appf2 = strcat(str_appf2, '_OBJ2_RUN3');
%str_appf2 = strcat(str_appf2, num2str(i));
str_appf2 = strcat(str_appf2, '.txt');
UML = importdata(str_appf2);
%figure(10+floor(rand*1000));
for j=1:20
%plot(a,UML(j,:),'-');
%x=(i-20000)/5000+1;
s=size(UML);
si=s(1,2);
b(i,j)=UML(j,si);
%h(j)=max(UML(j,:));
%hold on;
end
%title(test_case);xlabel('Time');ylabel('Upstreamwaterlevel');
end

% str = 'E:\pMOEAD for RFC\perfectPareto\';
% str2 = strcat(str, 'PF_');
% str2 = strcat(str2, test_case);
% str2 = strcat(str2, '.txt');
% PF = importdata(str2);
% plot(PF(:,1),PF(:,2),'b.')
% hold on


%  PF = importdata('E:\pMOEAD for RFCfu\perfectPareto\PF_ankang20001012.txt');
%  plot(PF(:,1),PF(:,2),'b.')
%  hold on


str_appf = 'E:\pMOEAD for RFCfu\POFmid\';

%str_appf = 'F:\xr\1186285856\FileRecv\code\pMOEAD for RFC\POF\';
total=0;
%for i=21000:1500:49500
for i=84000:6000:198000    

    %figure(i);

% str = 'E:\pMOEAD for RFC\perfectPareto\';
% str2 = strcat(str, 'PF_');
% str2 = strcat(str2, test_case);
% str2 = strcat(str2, '.txt');
% PF = importdata(str2);
PF = importdata('E:\pMOEAD for RFCfu\perfectPareto\PF_ankang20001012.txt');
plot(PF(:,1),PF(:,2),'b.')


% PF = importdata('E:\pMOEAD for RFCfu\perfectPareto\PF_ankang20100715.txt');
% plot(PF(:,1),PF(:,2),'b.')

s = strcat(num2str(i),test_case);
title(s);xlabel('Maximum Upstream Water Level (m)');ylabel('Maximum Out Flow (m3/s)');
hold on

str_appf2 = strcat(str_appf, num2str(i));
str_appf2 = strcat(str_appf2, '_POF_MOEAD_');
str_appf2 = strcat(str_appf2, test_case);
str_appf2 = strcat(str_appf2, '_OBJ2_RUN1');
%str_appf2 = strcat(str_appf2, num2str(i));
str_appf2 = strcat(str_appf2, '.txt');
POF = importdata(str_appf2);
%figure(10+floor(rand*1000));
plot(POF(:,1),POF(:,2),'gO')
%hold on 
  % x=(i-20000)/5000+1;
   M=[];num=0;
    for j=1:20
        if b(i,j)>=324 && b(i,j)<=326
           num=num+1;
           M(num,:)= POF(j,:); 
           plot(M(:,1),M(:,2),'rH')
        end    
    end 
    total=total+num
    if num>1
        maxa=max(M);
        MAX(i,:)=maxa;
        mina=min(M);
        MIN(i,:)=mina;
    elseif num==1
         MAX(i,:)=M;
         MIN(i,:)=M;
    end
    %axis([310 330 0 40000]);
      pause(2);
     hold off;
end
%title(test_case);xlabel('Maximum Upstream Water Level (m)');ylabel('Maximum Out Flow (m3/s)');
meanmaxpoint(1,1)=mean(nonzeros(MAX(:,1)))
meanmaxpoint(1,2)=mean(nonzeros(MAX(:,2)))
meanminpoint(1,1)=mean(nonzeros(MIN(:,1)))
meanminpoint(1,2)=mean(nonzeros(MIN(:,2)))
XR=(MAX+MIN)/2;
mid=(meanmaxpoint+meanminpoint)/2
