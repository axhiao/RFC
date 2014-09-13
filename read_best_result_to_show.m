function read_best_result_to_show
clc;
str_pr = 'PF\';
str_appf = 'M:\';
test_case= 'UF8';
index_best_result = 28; % the index of result with the minimum IGD.
num_prefer =2;        % the number of prefer areas.
%% PF
str_pr2 = strcat(str_pr, test_case);
str_pr2 = strcat(str_pr2, '.dat');
pftrue = importdata(str_pr2);
%% approximate solutions
str_appf2 = strcat(str_appf, 'POF_MOEAD_');
str_appf2 = strcat(str_appf2, test_case);
str_appf2 = strcat(str_appf2, '_RUN');
str_appf2 = strcat(str_appf2, num2str(index_best_result));
str_appf2 = strcat(str_appf2, '.txt');
Appro_pftrue = importdata(str_appf2);
%% plot the outline of problems
plot_outline(test_case);
%% show the result
%figure(10+floor(rand*1000));
if num_prefer == 1
    if strcmp(test_case, 'ZDT1') |   strcmp(test_case, 'ZDT4') ...
         | strcmp(test_case, 'UF1') | strcmp(test_case, 'UF2')...
         |strcmp(test_case, 'UF3')  
        point=[[0 0];[1 1]];
        outrank1= [[0 0];[3.3742300e-001  4.1911877e-001]];
        outrank2= [[0 0];[4.2121500e-001  3.5098921e-001]];
        plot(pftrue(:,1),pftrue(:,2),'k-',Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.',...
            point(:,1),point(:,2),'g-',outrank1(:,1),outrank1(:,2),'m-',outrank2(:,1),outrank2(:,2),'m-');title(test_case,'fontsize',24);xlabel('f1','fontsize',24);ylabel('f2','fontsize',24);set(gca,'fontsize',24);
    elseif strcmp(test_case, 'ZDT2')  |strcmp(test_case, 'ZDT6') |strcmp(test_case, 'UF4') 
        point=[[0 0];[1 1]];
        outrank1= [[0 0];[5.860880000000001e-001    6.565008562560000e-001]];
        outrank2= [[0 0];[6.490110000000000e-001    5.787847218790000e-001]];
        plot(pftrue(:,1),pftrue(:,2),'k-',Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.',...
            point(:,1),point(:,2),'g-',outrank1(:,1),outrank1(:,2),'c-',outrank2(:,1),outrank2(:,2),'c-');title(test_case,'fontsize',20);xlabel('f1','fontsize',20);ylabel('f2','fontsize',20);
    elseif strcmp(test_case, 'UF7')
         point=[[0 0];[1 1]];
        outrank1= [[0 0];[4.646446610000000e-001    5.353553390000000e-001]];
        outrank2= [[0 0];[5.353553390000000e-001    4.646446610000000e-001]];
         plot(pftrue(1:333,1),pftrue(1:333,2),'-k',pftrue(334:666,1),pftrue(334:666,2),'-k',pftrue(667:1000,1),pftrue(667:1000,2),'-k',Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.',...
            point(:,1),point(:,2),'g-',outrank1(:,1),outrank1(:,2),'c-',outrank2(:,1),outrank2(:,2),'c-');title(test_case,'fontsize',20);xlabel('f1','fontsize',20);ylabel('f2','fontsize',20);
    elseif strcmp(test_case, 'ZDT3')
        point=[[0 0];[1 1]];
        outrank1= [[0 0];[2.385465470000000e-001    2.883177207357462e-001]];
        outrank2= [[0 0];[2.577706300000000e-001    2.421610944958118e-001]];
        plot(pftrue(1:156,1),pftrue(1:156,2),'-k',pftrue(157:297,1),pftrue(157:297,2),'-k',pftrue(298:381,1),pftrue(298:381,2),'-k',pftrue(382:446,1),pftrue(382:446,2),'-k',pftrue(447:500,1),pftrue(447:500,2),'-k',...
            Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.',...
            point(:,1),point(:,2),'g-',outrank1(:,1),outrank1(:,2),'c-',outrank2(:,1),outrank2(:,2),'c-');title(test_case,'fontsize',20);xlabel('f1','fontsize',20);ylabel('f2','fontsize',20);
    elseif  strcmp(test_case, 'UF5')
        point=[[0 0];[1 1]];
        plot(pftrue(:,1),pftrue(:,2),'k.',Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.',...
            point(:,1),point(:,2),'g-');title(test_case,'fontsize',20);xlabel('f1','fontsize',20);ylabel('f2','fontsize',20);
    elseif strcmp(test_case, 'UF6')
        point=[[0 0];[0.6 1]];
        outrank1= [[0 0];[3.3964466e-001  6.6035534e-001]];
        outrank2= [[0 0];[4.1035534e-001  5.8964466e-001]];
         plot(pftrue(1:333,1),pftrue(1:333,2),'-k',pftrue(334:666,1),pftrue(334:666,2),'-k',pftrue(667:1000,1),pftrue(667:1000,2),'-k',Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.',...
            point(:,1),point(:,2),'g-',outrank1(:,1),outrank1(:,2),'c-',outrank2(:,1),outrank2(:,2),'c-');title(test_case,'fontsize',20);xlabel('f1','fontsize',20);ylabel('f2','fontsize',20);
    elseif strcmp(test_case, 'DTLZ1')
        point=[[0 0 0];[0.25 0.25 0.25]];
        plot3(Appro_pftrue(4:size(Appro_pftrue,1),1),Appro_pftrue(4:size(Appro_pftrue,1),2),Appro_pftrue(4:size(Appro_pftrue,1),3),'b.',point(:,1),point(:,2),point(:,3),'g-');
        title(test_case);xlabel('f1');ylabel('f2');
    elseif strcmp(test_case, 'DTLZ2') | strcmp(test_case, 'DTLZ3')  | strcmp(test_case, 'DTLZ4') 
        point=[[0 0 0];[1 1 1]];
        plot3(Appro_pftrue(4:size(Appro_pftrue,1),1),Appro_pftrue(4:size(Appro_pftrue,1),2),Appro_pftrue(4:size(Appro_pftrue,1),3),'b.',point(:,1),point(:,2),point(:,3),'g-');
        title(test_case);xlabel('f1');ylabel('f2');
    elseif strcmp(test_case, 'DTLZ6')
        point=[[0 0 0];[0.75 0.75 5]];
        plot3(Appro_pftrue(4:size(Appro_pftrue,1),1),Appro_pftrue(4:size(Appro_pftrue,1),2),Appro_pftrue(4:size(Appro_pftrue,1),3),'b.',point(:,1),point(:,2),point(:,3),'g-');
        title(test_case);xlabel('f1');ylabel('f2');
    elseif strcmp(test_case, 'UF8')
        point=[[0 0 0];[0.2222 0.4444 1]];
        plot3(Appro_pftrue(4:size(Appro_pftrue,1),1),Appro_pftrue(4:size(Appro_pftrue,1),2),Appro_pftrue(4:size(Appro_pftrue,1),3),'b.',point(:,1),point(:,2),point(:,3),'g-');
        title(test_case);xlabel('f1');ylabel('f2');
    elseif strcmp(test_case, 'UF9')
        point=[[0 0 0];[0.65 0.1 0.4]];
        plot3(Appro_pftrue(4:size(Appro_pftrue,1),1),Appro_pftrue(4:size(Appro_pftrue,1),2),Appro_pftrue(4:size(Appro_pftrue,1),3),'b.',point(:,1),point(:,2),point(:,3),'g-');
        title(test_case);xlabel('f1');ylabel('f2');
    elseif strcmp(test_case, 'UF10')
        point=[[0 0 0];[0.5  0.11111  1]];
        plot3(Appro_pftrue(4:size(Appro_pftrue,1),1),Appro_pftrue(4:size(Appro_pftrue,1),2),Appro_pftrue(4:size(Appro_pftrue,1),3),'b.',point(:,1),point(:,2),point(:,3),'g-');
        title(test_case);xlabel('f1');ylabel('f2');
    end
elseif num_prefer == 2
    if strcmp(test_case, 'ZDT1') | strcmp(test_case, 'ZDT4')|strcmp(test_case, 'UF1') ...
        | strcmp(test_case, 'UF2')  |strcmp(test_case, 'UF3')
        point=[[0 0];[1 0.5]];
        point2=[[0 0];[0.5 1]];
        outrank1= [[0 0];[2.152985008952569e-001    5.359973050775924e-001]];
        outrank2= [[0 0];[2.859538457824880e-001    4.652534752029818e-001]];
        outrank3= [[0 0];[4.948715695758288e-001    2.965289134755936e-001]];
        outrank4= [[0 0];[5.774317323384965e-001    2.401107104725738e-001]];
        plot(pftrue(:,1),pftrue(:,2),'k-',Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.',...
            point(:,1),point(:,2),'g-',point2(:,1),point2(:,2),'g-',outrank1(:,1),outrank1(:,2),'c-',outrank2(:,1),outrank2(:,2),'c-',...
            outrank3(:,1),outrank3(:,2),'c-',outrank4(:,1),outrank4(:,2),'c-');
        title(test_case);xlabel('f1');ylabel('f2');
    elseif strcmp(test_case, 'ZDT2')  |  strcmp(test_case, 'ZDT6') | strcmp(test_case, 'UF4') 
        point=[[0 0];[1 0.5]];
        point2=[[0 0];[0.5 1]];
        outrank1= [[0 0];[3.749637674949956e-001    8.594021730659589e-001]];
        outrank2= [[0 0];[4.520063713392312e-001    7.956902402687410e-001]];
        outrank3= [[0 0];[7.534744805314197e-001    4.322762071879072e-001]];
        outrank4= [[0 0];[8.074175664650349e-001    3.480768733636809e-001]];
        plot(pftrue(:,1),pftrue(:,2),'k-',Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.',...
            point(:,1),point(:,2),'g-',point2(:,1),point2(:,2),'g-',outrank1(:,1),outrank1(:,2),'c-',outrank2(:,1),outrank2(:,2),'c-',...
            outrank3(:,1),outrank3(:,2),'c-',outrank4(:,1),outrank4(:,2),'c-');
        title(test_case);xlabel('f1');ylabel('f2');    
    elseif  strcmp(test_case, 'ZDT3')
        point=[[0 0];[1 0.5]];
        point2=[[0 0];[0.5 1]];
        outrank1= [[0 0];[2.094067324424454e-001    4.814028914385481e-001]];
        outrank2= [[0 0];[2.227146999425844e-001    3.822927989209846e-001]];
        outrank3= [[0 0];[4.094319400000000e-001    2.405788180526862e-001]];
        outrank4= [[0 0];[4.158777369722964e-001    1.561648168587742e-001]];
        plot(pftrue(1:156,1),pftrue(1:156,2),'-k',pftrue(157:297,1),pftrue(157:297,2),'-k',pftrue(298:381,1),pftrue(298:381,2),'-k',pftrue(382:446,1),pftrue(382:446,2),'-k',pftrue(447:500,1),pftrue(447:500,2),'-k',...
            Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.',point(:,1),point(:,2),'g-',point2(:,1),point2(:,2),'g-',...
            outrank1(:,1),outrank1(:,2),'c-',outrank2(:,1),outrank2(:,2),'c-',...
            outrank3(:,1),outrank3(:,2),'c-',outrank4(:,1),outrank4(:,2),'c-');
        title(test_case);xlabel('f1');ylabel('f2');
    elseif  strcmp(test_case, 'UF5') 
        point=[[0 0];[0.42857 1]];
        point2=[[0 0];[1 0.42857]];        
        plot(pftrue(:,1),pftrue(:,2),'k.',Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.',...
            point(:,1),point(:,2),'g-',point2(:,1),point2(:,2),'g-');
    elseif strcmp(test_case, 'UF6')
        point=[[0 0];[0.6 1]];
        point2=[[0 0];[1 0.1428]];
        outrank1= [[0 0];[3.396446600000000e-001    6.603553400000000e-001]];
        outrank2= [[0 0];[4.103553390000000e-001    5.896446610000000e-001]];
        outrank3= [[0 0];[8.396884130413690e-001    1.603115869586310e-001]];
        outrank4= [[0 0];[9.103990913338491e-001    8.960090866615089e-002]];        
        plot(pftrue(1:333,1),pftrue(1:333,2),'-k',pftrue(334:666,1),pftrue(334:666,2),'-k',pftrue(667:1000,1),pftrue(667:1000,2),'-k',Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.',...
            point(:,1),point(:,2),'g-',point2(:,1),point2(:,2),'g-', outrank1(:,1),outrank1(:,2),'c-',outrank2(:,1),outrank2(:,2),'c-',...
            outrank3(:,1),outrank3(:,2),'c-',outrank4(:,1),outrank4(:,2),'c-');
        title(test_case);xlabel('f1');ylabel('f2');
    elseif strcmp(test_case, 'UF7')
        point=[[0 0];[1 0.5]];
        point2=[[0 0];[0.5 1]];
        outrank1= [[0 0];[2.979779941870930e-001    7.020220058129070e-001]];
        outrank2= [[0 0];[3.686886724795737e-001    6.313113275204263e-001]];
        outrank3= [[0 0];[6.313113275204263e-001    3.686886724795737e-001]];
        outrank4= [[0 0];[7.020220058129070e-001    2.979779941870930e-001]];
        plot(pftrue(:,1),pftrue(:,2),'k-',Appro_pftrue(3:size(Appro_pftrue,1),1),Appro_pftrue(3:size(Appro_pftrue,1),2),'b.',...
            point(:,1),point(:,2),'g-',point2(:,1),point2(:,2),'g-',outrank1(:,1),outrank1(:,2),'c-',outrank2(:,1),outrank2(:,2),'c-',...
            outrank3(:,1),outrank3(:,2),'c-',outrank4(:,1),outrank4(:,2),'c-');
        title(test_case);xlabel('f1');ylabel('f2');
    elseif strcmp(test_case, 'DTLZ1')
        point1=[[0 0 0];[0.5 0.1 0.3]];
        point2 =[[0 0 0];[0.3 0.5 0.1]];
        plot3(Appro_pftrue(4:size(Appro_pftrue,1),1),Appro_pftrue(4:size(Appro_pftrue,1),2),Appro_pftrue(4:size(Appro_pftrue,1),3),'b.',point1(:,1),point1(:,2),point1(:,3),'g-',point2(:,1),point2(:,2),point2(:,3),'g-');
    title(test_case);xlabel('f1');ylabel('f2');zlabel('f3');
    elseif strcmp(test_case, 'DTLZ2') | strcmp(test_case, 'DTLZ3')|strcmp(test_case, 'DTLZ4')
        point1=[[0 0 0];[1 0.2 0.6]];
        point2 =[[0 0 0];[0.6 1 0.2]];
        plot3(Appro_pftrue(4:size(Appro_pftrue,1),1),Appro_pftrue(4:size(Appro_pftrue,1),2),Appro_pftrue(4:size(Appro_pftrue,1),3),'b.',point1(:,1),point1(:,2),point1(:,3),'g-',point2(:,1),point2(:,2),point2(:,3),'g-');
    title(test_case);xlabel('f1');ylabel('f2');zlabel('f3');
    elseif strcmp(test_case, 'DTLZ6')
        point1=[[0 0 0];[0.75 0.75 5]];
        point2 =[[0 0 0];[0.75 0.13 5.3]];
        plot3(Appro_pftrue(4:size(Appro_pftrue,1),1),Appro_pftrue(4:size(Appro_pftrue,1),2),Appro_pftrue(4:size(Appro_pftrue,1),3),'b.',point1(:,1),point1(:,2),point1(:,3),'g-',point2(:,1),point2(:,2),point2(:,3),'g-');
    title(test_case);xlabel('f1');ylabel('f2');zlabel('f3');
    elseif strcmp(test_case, 'UF8')
        point1=[[0 0 0];[1 0.2 0.6]];
        point2 =[[0 0 0];[0.6 1 0.2]];
        plot3(Appro_pftrue(4:size(Appro_pftrue,1),1),Appro_pftrue(4:size(Appro_pftrue,1),2),Appro_pftrue(4:size(Appro_pftrue,1),3),'b.',point1(:,1),point1(:,2),point1(:,3),'g-',point2(:,1),point2(:,2),point2(:,3),'g-');
    title(test_case);xlabel('f1');ylabel('f2');zlabel('f3');
    elseif strcmp(test_case, 'UF9')
        point1=[[0 0 0];[0.65 0.1 0.4]];
        point2 =[[0 0 0];[0.1 0.65 0.4]]; 
        plot3(Appro_pftrue(4:size(Appro_pftrue,1),1),Appro_pftrue(4:size(Appro_pftrue,1),2),Appro_pftrue(4:size(Appro_pftrue,1),3),'b.',point1(:,1),point1(:,2),point1(:,3),'g-',point2(:,1),point2(:,2),point2(:,3),'g-');
    title(test_case);xlabel('f1');ylabel('f2');zlabel('f3');
    elseif strcmp(test_case, 'UF10')
        point1=[[0 0 0];[1 0.2 0.6]];
        point2 =[[0 0 0];[0.6 1 0.2]];
        plot3(Appro_pftrue(4:size(Appro_pftrue,1),1),Appro_pftrue(4:size(Appro_pftrue,1),2),Appro_pftrue(4:size(Appro_pftrue,1),3),'b.',point1(:,1),point1(:,2),point1(:,3),'g-',point2(:,1),point2(:,2),point2(:,3),'g-');
    title(test_case);xlabel('f1');ylabel('f2');zlabel('f3');
    end
end
end

function plot_outline(test_case)
if strcmp(test_case, 'DTLZ2') | strcmp(test_case, 'DTLZ3')  | strcmp(test_case, 'DTLZ4') ...
        | strcmp(test_case, 'UF8') | strcmp(test_case, 'UF10') 
        a=50;
        u=0:pi/a:2*pi;
        v=0:pi/a:2*pi;
        [U,V]=meshgrid(u,v); 
        x=sin(U).*cos(V);
        y=sin(U).*sin(V);
        z=cos(U);
        c=0*z;
        colormap(gray);
        mesh(abs(x),abs(y),abs(z),c);
        hold on;
elseif strcmp(test_case, 'DTLZ1')
        Data1 = [];Data2 = [];Data3 = [[0.5 0 0];[0 0.5 0]];
        N = 40;
        for i=0:1/N:0.5
            Data1 = [[i 0 0.5-i];[i 0.5-i 0]];
            Data2 = [[0 0.5-i i];[i 0.5-i 0]];
            plot3(Data1(:,1),Data1(:,2),Data1(:,3),'k-',Data2(:,1),Data2(:,2),Data2(:,3),'k-');
            hold on;
        end
        plot3(Data3(:,1),Data3(:,2),Data3(:,3),'k-');
        hold on;
elseif strcmp(test_case, 'DTLZ6')
        [X,Y] = meshgrid([0:0.02:1]);
        F=[];
        for i=1:size(X,1)
            for j=1:size(X,2)
                F(i,j) = 2*(3-X(i,j)/2*(1+sin(3*pi*X(i,j)))-Y(i,j)/2*(1+sin(3*pi*Y(i,j))));
            end
        end
        colormap(gray);
        mesh(X,Y,F,0*F);
        hold on;
elseif strcmp(test_case, 'UF9')
        Data1 = [];Data2 = []
        for i=0:0.05:0.25
            Data1 = [[0 0 1];[1-i i 0]];
            Data2 = [[0 0 1];[i 1-i 0]];
            plot3(Data1(:,1),Data1(:,2),Data1(:,3),'k-',Data2(:,1),Data2(:,2),Data2(:,3),'k-');
            hold on;
        end
        for i=0:0.1:1
            Data1 = [[i 0 1-i];[0.75*i 0.25*i 1-i]];
            Data2 = [[0.25*i 0.75*i 1-i];[0 i 1-i]];
            plot3(Data1(:,1),Data1(:,2),Data1(:,3),'k-',Data2(:,1),Data2(:,2),Data2(:,3),'k-');
            hold on;
        end
end    
end

