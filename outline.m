clc;
%% DTLZ2 DTLZ3 DTLZ4 UF8 UF10
% a=50;
% u=0:pi/a:2*pi;
% v=0:pi/a:2*pi;
% [U,V]=meshgrid(u,v); 
% x=sin(U).*cos(V);
% y=sin(U).*sin(V);
% z=cos(U);
% c=0*z;
% colormap(gray);
% mesh(abs(x),abs(y),abs(z),c);
%% DTLZ1
% Data1 = [];Data2 = [];Data3 = [[0.5 0 0];[0 0.5 0]];
% N = 20;
% for i=0:1/N:0.5
%     Data1 = [[i 0 0.5-i];[i 0.5-i 0]];
%     Data2 = [[0 0.5-i i];[i 0.5-i 0]];
%     plot3(Data1(:,1),Data1(:,2),Data1(:,3),'k-',Data2(:,1),Data2(:,2),Data2(:,3),'k-');
%    % xlabel('f1');ylalbel('f2');zlabel('f3');
%     hold on;
% end
% plot3(Data3(:,1),Data3(:,2),Data3(:,3),'k-');
%% DTLZ6
% [X,Y] = meshgrid([0:0.02:1]);
% F=[];
% for i=1:size(X,1)
%     for j=1:size(X,2)
%         F(i,j) = 2*(3-X(i,j)/2*(1+sin(3*pi*X(i,j)))-Y(i,j)/2*(1+sin(3*pi*Y(i,j))));
%     end
% end
% colormap(gray);
% mesh(X,Y,F,0*F);
% hold on;
%% UF9
% Data1 = [];Data2 = [];Data3 = [[0.5 0 0];[0 0.5 0]];
% for i=0:0.05:0.25
%     Data1 = [[0 0 1];[1-i i 0]];
%     Data2 = [[0 0 1];[i 1-i 0]];
%     plot3(Data1(:,1),Data1(:,2),Data1(:,3),'k-',Data2(:,1),Data2(:,2),Data2(:,3),'k-');
%     hold on;
% end
% for i=0:0.1:1
%     Data1 = [[i 0 1-i];[0.75*i 0.25*i 1-i]];
%     Data2 = [[0.25*i 0.75*i 1-i];[0 i 1-i]];
%     plot3(Data1(:,1),Data1(:,2),Data1(:,3),'k-',Data2(:,1),Data2(:,2),Data2(:,3),'k-');
%     hold on;
% end
% % str_pr = 'PF\';
% % test_case= 'UF9';
% % str_pr = strcat(str_pr, test_case);
% % str_pr = strcat(str_pr, '.dat');
% % pftrue = importdata(str_pr);
% % plot3(pftrue(:,1),pftrue(:,2),pftrue(:,3),'b.');

