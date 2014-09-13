function WS
str_pf = 'ParameterSetting\Weight vectr(MUD)\';
str_pf_new = 'ParameterSetting\Weight vector (WS(MUD))\';
m=15;
pop = 600;
str_pf2 = strcat('W', num2str(m));
str_pf2 = strcat(str_pf2, '_');
str_pf2 = strcat(str_pf2, num2str(pop));
test_instance = strcat(str_pf2, '.dat');

str_pf2 = strcat(str_pf, test_instance);
str_pf_new = strcat(str_pf_new, test_instance);
weight = importdata(str_pf2);
weight_new =[];
emxitong = 1e-6;
for i=1:pop
    multi = 1.0;
    for j=1:m
        multi = multi * weight(i,j);
    end
    if multi == 0
        weight(i,:) = weight(i,:) + emxitong;
    end
    weight_new = [weight_new; 1./weight(i,:)/(sum(1./weight(i,:)))];
end
save(str_pf_new,'weight_new','-ASCII');
