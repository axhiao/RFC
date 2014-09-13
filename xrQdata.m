str_appf = 'Qdata\';
%test_case='ankang20000712';
%test_case='ankang20030828';
%test_case='ankang20051001';
test_case='ankang20100715';

str_appf2 = strcat(str_appf, test_case);
str_appf2 = strcat(str_appf2, '.txt');
X = importdata(str_appf2);
a=0:1:max(size(X))-1;
plot(a,X(:,1))
title(test_case);xlabel('Time');ylabel('Input flow');