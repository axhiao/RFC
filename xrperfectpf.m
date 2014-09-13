str_appf = 'perfectPareto\';
%test_case='ankang20000712';
test_case='ankang20030828';
%test_case='ankang20051001';
%test_case='ankang20100715';
str_appf2 = strcat(str_appf, 'PF_');
str_appf2 = strcat(str_appf2, test_case);
str_appf2 = strcat(str_appf2, '.txt');
PF = importdata(str_appf2);
plot(PF(:,1),PF(:,2),'.')
hold on
