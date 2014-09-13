str_pf = 'PF\';
str_appf = 'POF\';
test_case= 'DTLZ2';
% A = [10 60];
% B = [25 20];
% C = [75 20];
% D = [85 60];
% E = [50 90];
% F = [70 80];
% point = [A;B;C;D;E;A];

% pftrue = importdata(str_pf);
% plot3(pftrue(:,10),pftrue(:,11),pftrue(:,12),'r.');
% % count = 5;
% % 
str_pf = strcat(str_pf, test_case);
str_pf = strcat(str_pf, '.dat');
pftrue = importdata(str_pf);
% 
figure(1);
for i=0:10:300
    str_pf2 = strcat(str_appf, num2str(i));
    str_pf2 = strcat(str_pf2, '\');
    str_pf2 = strcat(str_pf2, 'POF_MOEAD_');
    str_pf2 = strcat(str_pf2, test_case);
    str_pf2 = strcat(str_pf2, '_RUN');
    str_pf2 = strcat(str_pf2, num2str(2));
    str_pf2 = strcat(str_pf2, '.txt');
    appftrue = importdata(str_pf2);
    %figure(i);
    plot3(pftrue(:,1),pftrue(:,2),pftrue(:,3),'k.',appftrue(:,1),appftrue(:,2),appftrue(:,3),'b.');xlabel('f1', 'FontSize',20');    ylabel('f2', 'FontSize',20); zlabel('f3','FontSize',20);
    %plot(pftrue(:,1),pftrue(:,2),'k-',appftrue(:,1),appftrue(:,2),'b.');  % axis([0 100 0 100]); set(gca,'fontsize',20);  xlabel('X1', 'FontSize',20');    ylabel('X2', 'FontSize',20);
    %h=legend('MOEA/D');
    %h=legend('NIW-MOEA/D');

%    axis([0 1 0 5]);
%     plot3(pftrue(:,1),pftrue(:,2),pftrue(:,3),'r.');
%     axis([0 1.2 0 1.2 0 1.2]);
    %xlabel('f1');
    pause(0.01);
end
% pf = 'F:\maxiaoliang\≥£”√\PFture\';
% pf = strcat(pf, test_case);
% pf = strcat(pf, '.dat');
% p = importdata(pf);


% plot3(pftrue(:,1),pftrue(:,2),pftrue(:,3),'r.');
%plot(pftrue(:,1),pftrue(:,2),'r.');
%plot3(pftrue(:,1),pftrue(:,2),pftrue(:,3),'r*',p(:,1),p(:,2),p(:,3),'k.');