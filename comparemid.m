clc;clear all;close all;
%16\14\18\24
%T=6,3,4,6
 clc;clear all;
 a=6:6:6*16;
 %a=3:3:3*14;
 %a=4:4:4*18;
 %a=6:6:6*24;

str_appf = 'E:\421\2ndmid516\UpstreamWaterLevelmid\';
%str_appf = 'E:\421\pMOEAD for RFC\UpstreamWaterLevelmid\';
test_case='ankang20000712';
%test_case='ankang20030828';
%test_case='ankang20051001';
%test_case='ankang20100715';
%count = 50000;
for  k=1:10
for i=20000:5000:50000
str_appf2 = strcat(str_appf, num2str(i));
str_appf2 = strcat(str_appf2, '_UWL_MOEAD_');
str_appf2 = strcat(str_appf2, test_case);
str_appf2 = strcat(str_appf2, '_OBJ2_RUN');
str_appf2 = strcat(str_appf2, num2str(k));
str_appf2 = strcat(str_appf2, '.txt');
UML = importdata(str_appf2);
%figure(10+floor(rand*1000));
for j=1:20
plot(a,UML(j,:),'-');
x=(i-20000)/5000+1;
s=size(UML);
si=s(1,2);
b(x,j)=UML(j,si);
%h(j)=max(UML(j,:));
hold on;
end
hold on;
end
end
title(test_case);xlabel('Time');ylabel('Upstreamwaterlevel');
figure

% str_appf = 'E:\421\pMOEAD for RFC\515\perfectPareto\';
% str_appf2 = strcat(str_appf, 'PF_');
% str_appf2 = strcat(str_appf2, test_case);
% str_appf2 = strcat(str_appf2, '.txt');
% PF = importdata(str_appf2);
% plot(PF(:,1),PF(:,2),'.')
% hold on

str_appf = 'E:\421\2ndmid516\POFmid\';
%str_appf = 'E:\421\pMOEAD for RFC\POFmid\';
for k=1:10
total=0;
for i=20000:5000:50000;
str_appf2 = strcat(str_appf, num2str(i));
str_appf2 = strcat(str_appf2, '_POF_MOEAD_');
str_appf2 = strcat(str_appf2, test_case);
str_appf2 = strcat(str_appf2, '_OBJ2_RUN');
str_appf2 = strcat(str_appf2, num2str(k));
str_appf2 = strcat(str_appf2, '.txt');
POF = importdata(str_appf2);
%figure(10+floor(rand*1000));
plot(POF(:,1),POF(:,2),'yO')
%hold on 
   x=(i-20000)/5000+1;
   M=[];num=0;
    for j=1:20
        if b(x,j)>324 && b(x,j)<326
           num=num+1;
           M(num,:)= POF(j,:); 
           plot(M(:,1),M(:,2),'rH')
        end    
    end 
    total=total+num;
    kl(k)=total;
    if num>1
        maxa=max(M);
        MAX(x,:)=maxa;
        mina=min(M);
        MIN(x,:)=mina;
    elseif num==1
         MAX(x,:)=M;
         MIN(x,:)=M;
    end
    
  % figure 
end
end
title(test_case);xlabel('Maximum Upstream Water Level (m)');ylabel('Maximum Out Flow (m3/s)');
meanmaxpoint(1,1)=mean(nonzeros(MAX(:,1)));
meanmaxpoint(1,2)=mean(nonzeros(MAX(:,2)));
meanminpoint(1,1)=mean(nonzeros(MIN(:,1)));
meanminpoint(1,2)=mean(nonzeros(MIN(:,2)));
% if num>1
%         max=max(M);
%     end

% for i=1:(50000-20000)/5000+1
%     M=[];num=0;
%     for j=1:20
%         if b(i,j)>324 && b(i,j)<326
%            num=num+1;
%            M(num,:)= POF(j,:); 
%            plot(M(:,1),M(:,2),'rH')
%         end        
%     end 
%     
% end

% if num==0 
%     break   
% else
%   sizem=size(M);
%   l=sizem(1,1);w=sizem(1,2);
%   max=max(M);
%   min=min(M);
%   x0=0.5*(max(1,1)+min(1,1));
%   y0=0.5*(max(1,2)+min(1,2));
%   R=0.5*norm(max-min);
%   %circle(R,x0,y0);
%    r=1.1*R ;
%   %axis([310 330 0 40000])
%    %set(gca,'clipping','on')
% end


% if num==1
%         max=M;
%         min=M;
%     elseif num>1
%          max=max(M);
%          min=min(M);
%          x0(i)=0.5*(max(1,1)+min(1,1));
%          y0(i)=0.5*(max(1,2)+min(1,2));
%          v1(i)=0.5*(max(1,1)-min(1,1));
%          v2(i)=0.5*(max(1,2)-min(1,2));
% end

% if num==1
%         max=M;
%         min=M;
%     elseif num>1
%          max=max(M);
%          min=min(M);
%          x0(i)=0.5*(max(1,1)+min(1,1));
%          y0(i)=0.5*(max(1,2)+min(1,2));
%          v1(i)=0.5*(max(1,1)-min(1,1));
%          v2(i)=0.5*(max(1,2)-min(1,2));
%  end
%     




